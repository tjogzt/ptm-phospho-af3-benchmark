#!/usr/bin/env python3
"""
run_esmfold.py — 磷酸化感知 ESMFold 推理（B7）
======================================================================
Phosphorylation-aware ESMFold inference for the PTM benchmark (B7).

中文（Chinese）
---------------
目的：为论文 Fig 2 提供**第二条对比线**（ESMFold），检验"AF3 无响应"是否是
模型特异的还是普遍现象。

关键设计（保证与 AF3 线可比）：
  1. **序列同源**：ESMFold 直接使用 AF3 unphos 模型 CIF 中记录的全长 UniProt
     序列（残基编号 = UniProt 坐标），确保两条对比线在完全相同的序列上进行。
     （注意：`data/fasta/*.fasta` 为 v1 陈旧数据，切勿使用。）
  2. **磷酸化感知**：ESMFold **不解析 CCD_ 修饰残基**（SEP/TPO/PTR），无法像
     AF3 那样显式输入磷酸基团。故采用文献常用的 **phosphomimetic** 策略：
     在目标位点把 pSer/pThr → Asp(D)、pTyr → Glu(E)，以负电荷近似磷酸基团的
     静电与空间效应。
  3. **两态协议一致**：unphos = 野生型序列；phos = 位点 phosphomimetic。两态
     各推理一次（ESMFold 为确定性模型，无 seed×model 采样），保存 PDB 供
     下游 compare_esmfold_vs_pdb.py 与实验/ AF3 位移比较。

局限（须写入论文 §3.3 / §5）：phosphomimetic 只近似电荷，不含真实磷酸基团
的立体与氢键特征；ESMFold 亦无 MSA 输入。故此线只能回答"另一类单序列预测器
是否对位点突变产生方向性响应"，不能等同于真实磷酸化建模。

English
-------
Goal: provide the second comparison line (ESMFold) for Fig 2, testing whether
the AF3 non-response is model-specific or generic.
  (1) SEQUENCE-INHERITED: ESMFold uses the full-length UniProt sequence recorded
      in the AF3 unphos CIF (residue numbering = UniProt), so both lines run on
      identical sequences. (Do NOT use the stale v1 `data/fasta/*.fasta`.)
  (2) PHOSPHO-AWARE: ESMFold does not parse CCD_ modifications, so we apply the
      standard phosphomimetic substitution (pS/pT -> Asp, pY -> Glu).
  (3) CONSISTENT TWO-STATE PROTOCOL: unphos = wild type; phos = phosphomimetic.
      ESMFold is deterministic (no seed x model sampling); PDBs are saved for the
      downstream comparison script.

输出 / Output:
  results/esm/pdb/<gene>_<state>_esmfold.pdb
  results/esm/esmfold_summary.json

用法 / Usage:
  python3 scripts/esm/run_esmfold.py [--device mps|cpu] [--limit N] [--force]
"""
import os
import sys
import json
import argparse
import time
from pathlib import Path

import numpy as np
import torch

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
sys.path.insert(0, str(_ROOT / "scripts" / "analysis"))
from benchmark_dataset import BENCHMARK_PROTEINS
from struct_utils import load_ca

CIFDIR = _ROOT / "results" / "af3" / "models" / "v2"
OUT = _ROOT / "results" / "esm"
PDBDIR = OUT / "pdb"

# ESMFold v1 权重来源：本机 huggingface.co 被 DNS 劫持，改用 ModelScope 官方镜像
# （facebook/esmfold_v1，同一份 HF 格式权重）。可用 ESMFOLD_PATH 覆盖。
# ESMFold v1 weights: huggingface.co is DNS-hijacked on this machine, so we use
# the ModelScope mirror of the identical HF checkpoint (override via ESMFOLD_PATH).
ESM_MODEL_PATH = os.environ.get("ESMFOLD_PATH",
                                "/Volumes/thinkplus/.model_cache/facebook/esmfold_v1")


def _model_ref():
    p = Path(ESM_MODEL_PATH)
    return str(p) if p.exists() else "facebook/esmfold_v1"

MANIFEST = _ROOT / "results" / "manifest_frozen.json"
_mf = json.loads(MANIFEST.read_text())
PRIMARY = list(_mf["cohorts"]["primary"])

# phosphomimetic 映射：磷酸化残基母体 → 替代残基
# Phosphomimetic map: parent residue -> substitution
PHOSPHOMIMETIC = {"S": "D", "T": "D", "Y": "E"}


def af3_sequence(gene, state):
    """从 AF3 模型 CIF 提取全长序列 + 残基编号（确保与 AF3 线同源）。

    Extract the full-length sequence and residue numbering from an AF3 model CIF
    so the ESMFold line runs on exactly the same sequence as the AF3 line.
    """
    d = CIFDIR / gene / state
    cand = sorted(p for p in d.glob("fold_*_model_0.cif") if "template" not in p.name)
    if not cand:
        return None, None
    _, seq, resids, _ = load_ca(cand[0], chain=None)
    return seq, list(resids)


def site_indices(seq, resids, sites):
    """把 UniProt 位点号映射为序列 0-based 索引。"""
    pos = {r: i for i, r in enumerate(resids)}
    return [pos[s] for s in sites if s in pos]


def make_phosphomimetic(seq, idxs):
    """在给定位点做 phosphomimetic 替换；返回新序列与替换记录。"""
    s = list(seq)
    changes = []
    for i in idxs:
        wt = s[i]
        mt = PHOSPHOMIMETIC.get(wt)
        if mt:
            s[i] = mt
            changes.append({"pos": i, "wt": wt, "mt": mt})
    return "".join(s), changes


def load_model(device):
    """加载 ESMFold（transformers 纯 PyTorch 实现，不依赖 openfold）。"""
    from transformers import AutoTokenizer, EsmForProteinFolding
    ref = _model_ref()
    print(f"Loading ESMFold v1 (transformers EsmForProteinFolding) from {ref}...")
    tok = AutoTokenizer.from_pretrained(ref)
    model = EsmForProteinFolding.from_pretrained(ref, low_cpu_mem_usage=True)
    # 注意：fp16 会使 transformers 的 compute_tm 崩溃——其在 loss.py 中用
    # `(weighted == max).nonzero()[0]`，低精度下无精确匹配 → 空张量 → IndexError。
    # 故统一使用 fp32（本机 128 GB 内存，显存/统一内存充裕）。
    # NOTE: fp16 breaks transformers compute_tm, which does
    # `(weighted == max).nonzero()[0]`; in reduced precision there is no exact
    # match, yielding an empty tensor and IndexError. Use fp32 throughout.
    model = model.float().eval()
    model = model.to(device)
    return tok, model


def fold_one(tok, model, seq, device):
    """对单条序列推理，返回 (pdb_str, plddt_list, ptm)。"""
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs)
    pdb = model.output_to_pdb(out)[0]
    plddt = out["plddt"][0].float().mean(-1).cpu().numpy().tolist()   # 每残基
    # 注意：批次为 1 时 transformers 返回的 ptm 是 0 维标量张量（shape ()），
    # 直接 [0] 索引会报 "index 0 is out of bounds for dimension 0 with size 0"。
    # NOTE: with batch size 1, transformers returns ptm as a 0-dim scalar tensor
    # (shape ()), so indexing [0] raises IndexError; flatten first.
    ptm_t = out["ptm"].float().cpu().reshape(-1)
    ptm = float(ptm_t[0]) if ptm_t.numel() else float("nan")
    # 清理缓存，避免 MPS 内存碎片累积
    del out
    if device == "mps":
        torch.mps.empty_cache()
    return pdb, plddt, ptm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default=None, choices=["mps", "cpu"],
                    help="默认优先 mps，不可用则 cpu")
    ap.add_argument("--limit", type=int, default=0, help="仅处理前 N 个蛋白（调试）")
    ap.add_argument("--force", action="store_true", help="覆盖已存在的 PDB")
    args = ap.parse_args()

    PDBDIR.mkdir(parents=True, exist_ok=True)

    if args.device:
        device = args.device
    else:
        device = "mps" if torch.backends.mps.is_available() else "cpu"
    print(f"Device: {device}  | torch {torch.__version__}")

    tok, model = load_model(device)

    cohort = PRIMARY[: args.limit] if args.limit else PRIMARY
    summary, t0 = [], time.time()

    for key in cohort:
        info = BENCHMARK_PROTEINS.get(key)
        if not info:
            print(f"[skip] {key}: not in benchmark")
            continue
        gene = info["gene"].lower()
        sites = list(info.get("phos_sites_uniprot") or
                     [int("".join(c for c in s if c.isdigit())) for s in info.get("phos_sites", [])])

        seq, resids = af3_sequence(gene, "unphos")
        if seq is None:
            print(f"[skip] {gene}: no AF3 CIF")
            continue
        idxs = site_indices(seq, resids, sites)
        if not idxs:
            print(f"[warn] {gene}: sites {sites} not found in AF3 numbering; "
                  f"resid range {resids[0]}..{resids[-1]}")
        phos_seq, changes = make_phosphomimetic(seq, idxs)

        rec = {"protein": key, "gene": gene, "length": len(seq),
               "sites_uniprot": sites, "phosphomimetic_changes": changes,
               "n_sites_mutated": len(changes)}
        print(f"\n=== {key} ({gene}) len={len(seq)} sites={sites} "
              f"mutations={len(changes)} ===")

        for state, s in (("unphos", seq), ("phos", phos_seq)):
            outp = PDBDIR / f"{gene}_{state}_esmfold.pdb"
            if outp.exists() and not args.force:
                print(f"  [cached] {outp.name}")
                # 仍记录指标（从摘要复用）
                rec[f"{state}_cached"] = True
                continue
            ts = time.time()
            try:
                pdb, plddt, ptm = fold_one(tok, model, s, device)
            except Exception as e:
                print(f"  [ERROR] {state}: {type(e).__name__}: {e}")
                rec[f"{state}_error"] = f"{type(e).__name__}: {e}"
                continue
            outp.write_text(pdb)
            rec[f"{state}_mean_plddt"] = round(float(np.mean(plddt)), 2)
            rec[f"{state}_ptm"] = round(ptm, 4)
            rec[f"{state}_plddt"] = [round(x, 2) for x in plddt]
            print(f"  {state}: pTM={ptm:.3f} pLDDT={np.mean(plddt):.1f} "
                  f"({time.time()-ts:.1f}s) -> {outp.name}")

        summary.append(rec)

    (OUT / "esmfold_summary.json").write_text(
        json.dumps({"device": device, "n_proteins": len(summary),
                    "phosphomimetic_rule": "pS/pT->D, pY->E",
                    "caveat": ("ESMFold 不解析 CCD_ 修饰；phos 态用 phosphomimetic 近似，"
                               "仅近似电荷效应，非真实磷酸化建模。"),
                    "proteins": summary}, indent=2, ensure_ascii=False))
    print(f"\n[OK] {OUT}/esmfold_summary.json  ({len(summary)} proteins, "
          f"{time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
