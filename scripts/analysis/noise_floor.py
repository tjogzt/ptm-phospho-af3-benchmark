#!/usr/bin/env python3
"""
noise_floor.py — B5：同态噪声地板对照
======================================================================
Same-state noise floor: does the phosphorylation-induced experimental
displacement exceed the displacement you observe between two structures
that are in the SAME phosphorylation state?

中文（Chinese）
---------------
主结果（C3/C4）是在"清洗后队列"上观察到 AF3 预测位移与实验位移无相关。
审稿人必问：实验位移（约 0.25–2.51 Å）是否本身就在噪声量级内？若在，则
阴性结果不可解释。本脚本用两类**同态**（same-state）对照回答该问题：

  (1) 实验同态复本地板（主）：同一晶体中若存在目标蛋白的**多条同源链**
      （非对称单元复本，NCS），将每条复本叠合到目标链后，逐残基 Cα 位移
      即为"同一状态、同一实验"下的测量离散 → 实验噪声地板上界。
      该量不涉及磷酸化，纯粹刻画结构测定的可重复性。

  (2) 预测器同态地板（次）：每个态下 AF3 的 15 个 seed×model 预测，两两
      （以第 0 号模型为参考）叠合后的逐残基离散 → 预测器在同一输入下的
      内部不重复度。用于说明"模型侧"的噪声量级。

English
-------
The primary result (C3/C4) is a null correlation on a purity-cleaned cohort.
A reviewer will ask whether the experimental displacements themselves are
within noise. We answer with two SAME-STATE controls:
  (1) experimental replicate floor: within one crystal, extra copies of the
      target protein (non-crystallographic symmetry mates) are superposed
      onto the reference chain; the per-residue Cα deviation is an empirical
      same-state, same-experiment measurement floor (upper bound, since NCS
      mates may also differ functionally).
  (2) predictor reproducibility floor: across the 15 AF3 seed×model
      predictions of one state, the per-residue spread after superposition.

输出 / Output:
  results/analysis/noise_floor.json
  results/figures/FigS1_noise_floor.pdf / .png

用法 / Usage:
  python3 scripts/analysis/noise_floor.py
"""
import sys, json
from pathlib import Path
from collections import defaultdict
from difflib import SequenceMatcher
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
sys.path.insert(0, str(_HERE))
from benchmark_dataset import BENCHMARK_PROTEINS
from struct_utils import aa_from_resname, load_ca

PDBDIR = _ROOT / "data" / "pdb"
CIFDIR = _ROOT / "results" / "af3" / "models" / "v2"
OUT = _ROOT / "results" / "analysis"
FIG = _ROOT / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)

MANIFEST = _ROOT / "results" / "manifest_frozen.json"
_mf = json.loads(MANIFEST.read_text())
PRIMARY = list(_mf["cohorts"]["primary"])

SHORT = {
    "CDK2": "CDK2", "ERK2": "ERK2", "Aurora_A": "AURKA", "IRF3": "IRF3",
    "P38alpha": "p38", "CDK7": "CDK7", "c-Src_kinase": "Src",
    "AMPK_alpha2": "AMPK", "Glycogen_Phosphorylase": "PYGM", "CDK1": "CDK1",
    "GSK3B": "GSK3B", "PAK1": "PAK1", "AKT1": "AKT1", "TBK1": "TBK1",
    "IKBKB": "IKBb",
}

_pdb = PDBParser(QUIET=True)
_cif = MMCIFParser(QUIET=True)


def kabsch(P, Q):
    """叠合 P 到 Q：P@R + t ≈ Q（行向量约定）。errstate 屏蔽 Accelerate 伪警告。"""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        H = Pc.T @ Qc
        U, _, Vt = np.linalg.svd(H)
        d = 1.0 if float(np.linalg.det(U @ Vt)) >= 0 else -1.0
        R = U @ np.diag([1.0, 1.0, d]) @ Vt
        t = Q.mean(0) - P.mean(0) @ R
    return R, t


def align_pairs(a, b):
    sm = SequenceMatcher(None, a, b, autojunk=False)
    out = []
    for blk in sm.get_matching_blocks():
        for k in range(blk.size):
            out.append((blk.a + k, blk.b + k))
    return out


def chain_ca(residues):
    """从残基列表提取 (coords[n,3], seq[str], resids[list])；纳入 SEP/TPO/PTR。"""
    coords, seq, resids = [], [], []
    for res in residues:
        if "CA" not in res:
            continue
        if res.id[0][0] not in (" ", "H"):
            continue
        aa = aa_from_resname(res.get_resname())
        if aa is None:
            continue
        coords.append(res["CA"].get_coord())
        seq.append(aa)
        resids.append(res.id[1])
    return (np.asarray(coords, float), "".join(seq), resids)


def all_chains(path):
    """返回 {chain_id: (coords, seq, resids)}（按文件格式自动选择解析器）。"""
    path = Path(path)
    parser = _cif if path.suffix.lower() in (".cif", ".mmcif") else _pdb
    s = parser.get_structure("x", str(path))
    model = next(iter(s))
    out = {}
    for ch in model:
        c = chain_ca(list(ch))
        if len(c[1]) >= 5:
            out[ch.get_id()] = c
    return out


def identity(a, b):
    """序列一致性 + 覆盖度（以较短者为分母）。"""
    m = align_pairs(a, b)
    if not m:
        return 0.0, 0.0
    ident = sum(1 for i, j in m if a[i] == b[j]) / max(len(a), len(b))
    cov = len(m) / min(len(a), len(b))
    return ident, cov


def replicate_displacements(ref, others, min_ident=0.95, min_cov=0.90):
    """把每个同源复本叠合到参考链，返回逐残基 |Δ| 列表（跨复本汇总）。

    Superpose each homologous replicate onto the reference chain and return
    the pooled per-residue |Δ| (across replicates).
    """
    rc, rseq, rres = ref
    pool = []
    n_rep = 0
    for cid, (cc, cseq, cres) in others.items():
        ident, cov = identity(rseq, cseq)
        if ident < min_ident or cov < min_cov:
            continue                      # 非同源链（如复合物伴侣蛋白）跳过
        m = align_pairs(rseq, cseq)
        if len(m) < 20:
            continue
        i_ref = np.array([x[0] for x in m])
        i_rep = np.array([x[1] for x in m])
        P = cc[i_rep]                     # 复本
        Q = rc[i_ref]                     # 参考
        with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
            R, t = kabsch(P, Q)
            dist = np.linalg.norm((P @ R + t) - Q, axis=1)
        pool.extend(dist.tolist())
        n_rep += 1
    return pool, n_rep


def model_floor(gene, state, allowed_resids=None, ref_idx=0):
    """同一态下 15 个 AF3 模型之间的逐残基离散（以第 ref_idx 号模型为参考）。

    Per-residue spread across AF3 models of one state (reference = model 0).

    allowed_resids: 仅统计落在实验结构（PDB）覆盖范围内的残基。AF3 预测的是
    全长 UniProt 序列，其无序末端在不同 seed 间差异极大；若不加限制，预测地板
    会被末端柔性淹没，与在**有序核心**上计算的实验位移不可比。
    Restrict to residues covered by the experimental structure: AF3 predicts the
    full-length sequence whose disordered termini differ wildly across seeds, so
    an unrestricted floor is not comparable to the experimental (ordered-core)
    displacement.
    """
    d = CIFDIR / gene / state
    files = sorted(d.glob("fold_*_seed*_model_*.cif")) if d.exists() else []
    if len(files) < 2:
        return None, 0
    loaded = []
    for f in files:
        coords, seq, resids, _ = load_ca(f, chain=None)
        loaded.append((coords, seq, resids))
    base = loaded[ref_idx]
    pool = []
    for k in range(len(loaded)):
        if k == ref_idx:
            continue
        cur = loaded[k]
        m = align_pairs(base[1], cur[1])
        if allowed_resids is not None:
            m = [(i, j) for (i, j) in m if base[2][i] in allowed_resids]
        if len(m) < 20:
            continue
        ib = np.array([x[0] for x in m])
        ic = np.array([x[1] for x in m])
        with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
            R, t = kabsch(cur[0][ic], base[0][ib])
            pool.extend(np.linalg.norm((cur[0][ic] @ R + t) - base[0][ib], axis=1).tolist())
    if not pool:
        return None, 0
    return pool, len(loaded) - 1


def main():
    exp_floor = []          # 每蛋白：实验同态复本地板
    pred_floor = []         # 每蛋白：预测器同态地板
    print("=" * 96)
    print("  B5 同态噪声地板：同一状态下可观测到的 Cα 位移量级")
    print("=" * 96)
    print(f"  {'蛋白':<8} {'实验同态均值':>10} {'p95':>7} {'复本数':>6} "
          f"{'预测同态均值':>10} {'p95':>7} {'模型数':>6}")

    for key in PRIMARY:
        info = BENCHMARK_PROTEINS.get(key)
        if not info or not info.get("paired"):
            continue
        gene = info["gene"].lower()
        short = SHORT.get(key, info["gene"])
        pdb_chains = info.get("pdb_chains", {})

        exp_pool, n_rep_tot, pdb_resids = [], 0, set()
        for state in ("unphos", "phos"):
            pid = next((x for x in info["pdbs"].get(state, [])
                        if (PDBDIR / f"{x.upper()}.pdb").exists()), None)
            if not pid:
                continue
            ch = pdb_chains.get(pid, pdb_chains.get(pid.upper(), info.get("chain", "A")))
            chains = all_chains(PDBDIR / f"{pid.upper()}.pdb")
            if ch not in chains:
                continue
            pdb_resids.update(chains[ch][2])          # 实验结构覆盖的有序残基
            others = {c: v for c, v in chains.items() if c != ch}
            pool, n_rep = replicate_displacements(chains[ch], others)
            if pool:
                exp_pool.extend(pool)
                n_rep_tot += n_rep

        pred_pool, n_mod = [], 0
        for state in ("unphos", "phos"):
            p, nm = model_floor(gene, state, allowed_resids=pdb_resids or None)
            if p:
                pred_pool.extend(p)
                n_mod = max(n_mod, nm)

        exp_stat = (round(float(np.mean(exp_pool)), 3),
                    round(float(np.percentile(exp_pool, 95)), 3), n_rep_tot) if exp_pool else (None, None, 0)
        pred_stat = (round(float(np.mean(pred_pool)), 3),
                     round(float(np.percentile(pred_pool, 95)), 3), n_mod) if pred_pool else (None, None, 0)

        print(f"  {short:<8} {str(exp_stat[0]):>10} {str(exp_stat[1]):>7} {exp_stat[2]:>6} "
              f"{str(pred_stat[0]):>10} {str(pred_stat[1]):>7} {pred_stat[2]:>6}")

        exp_floor.append({"protein": key, "short": short,
                          "exp_mean": exp_stat[0], "exp_p95": exp_stat[1], "n_replicates": exp_stat[2],
                          "pred_mean": pred_stat[0], "pred_p95": pred_stat[1], "n_models": pred_stat[2]})

    # ---- 与实验磷酸化位移对比（来自 compare 产物）----
    comp = json.loads((OUT / "af3_vs_pdb_summary.json").read_text())
    disp = {r["protein"]: r["mean_disp_exp"] for r in comp["per_protein"]
            if r.get("mean_disp_exp") is not None and not r.get("same_structure")}
    for row in exp_floor:
        row["phos_disp_exp"] = disp.get(row["protein"])

    exp_vals = [r["exp_mean"] for r in exp_floor if r["exp_mean"] is not None]
    pred_vals = [r["pred_mean"] for r in exp_floor if r["pred_mean"] is not None]
    disp_vals = [v for v in disp.values() if v is not None]

    summary = {
        "description": ("同一状态(same-state)下可观测位移量级；实验侧取同晶体同源复本(NCS)，"
                        "预测侧取同态 15 个 AF3 模型间的离散。用于判断实验磷酸化位移是否超出噪声。"),
        "n_proteins": len(exp_floor),
        "experimental_same_state_floor": {
            "median_mean": round(float(np.median(exp_vals)), 3) if exp_vals else None,
            "max_mean": round(float(np.max(exp_vals)), 3) if exp_vals else None,
            "n_proteins_with_replicates": len(exp_vals),
        },
        "predictor_same_state_floor": {
            "median_mean": round(float(np.median(pred_vals)), 3) if pred_vals else None,
            "max_mean": round(float(np.max(pred_vals)), 3) if pred_vals else None,
        },
        "phos_induced_experimental_displacement": {
            "median": round(float(np.median(disp_vals)), 3) if disp_vals else None,
            "range": [round(float(np.min(disp_vals)), 2), round(float(np.max(disp_vals)), 2)] if disp_vals else None,
        },
        "per_protein": exp_floor,
        "caveat": ("NCS 复本可能含真实构象异质性，故实验地板为噪声的**上界**（保守估计）；"
                   "预测地板反映采样不可重复度，非生物学变异。"),
    }
    (OUT / "noise_floor.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    # ---- 图（FigS1）----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 7.5,
                         "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "legend.fontsize": 6.3, "axes.linewidth": 0.6, "pdf.fonttype": 42})
    rows = [r for r in exp_floor if r["phos_disp_exp"] is not None and r["exp_mean"] is not None]
    rows.sort(key=lambda r: r["phos_disp_exp"])
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    y = np.arange(len(rows))
    ax.barh(y, [r["phos_disp_exp"] for r in rows], color="#4c72b0", height=0.55,
            label="Phospho-induced displacement", zorder=2)
    ax.barh(y, [r["exp_mean"] for r in rows], color="#dd8452", height=0.28,
            label="Same-state replicate floor", zorder=3)
    ax.set_yticks(y); ax.set_yticklabels([r["short"] for r in rows])
    ax.set_xlabel(r"Mean C$\alpha$ displacement ($\AA$)")
    ax.set_title("Same-state noise floor vs phospho-induced displacement", loc="left")
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(pad=0.6)
    fig.savefig(FIG / "FigS1_noise_floor.pdf", bbox_inches="tight")
    fig.savefig(FIG / "FigS1_noise_floor.png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    print("-" * 96)
    es = summary["experimental_same_state_floor"]
    ps = summary["predictor_same_state_floor"]
    ds = summary["phos_induced_experimental_displacement"]
    print(f"  实验同态地板（中位/最大，n={es['n_proteins_with_replicates']}）："
          f"{es['median_mean']} / {es['max_mean']} Å")
    print(f"  预测同态地板（中位/最大）：{ps['median_mean']} / {ps['max_mean']} Å")
    print(f"  实验磷酸化位移（中位 / 范围）：{ds['median']} Å / {ds['range']}")
    print("-" * 96)
    print(f"[OK] {OUT}/noise_floor.json")
    print(f"[OK] {FIG}/FigS1_noise_floor.pdf")
    return 0


if __name__ == "__main__":
    main()
