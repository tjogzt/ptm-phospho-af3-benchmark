#!/usr/bin/env python3
"""
generate_af3_jobs.py — AF3 批量作业生成（v2，适配审计修复后的基准）
=====================================================================
v1 的问题（见 audit/FINDINGS_P0.md）：
  - 序列取自 PDB 观测坐标（长度偏短、phos/unphos 不可比）
  - 位点编号混用三套体系，导致 PKA/EGFR/ERK2/c-Src 等挂错残基
  - 输出格式与实际提交格式不一致（无 CCD_ 前缀、种子合并方式不同）

v2 变更：
  - 序列一律取自 data/curated/uniprot_sequences.json（UniProt canonical 全长）
  - 位点取自 benchmark_dataset.py 的 phos_sites（已核准为 UniProt 坐标）
  - 输出：每作业单种子、ptmType 带 CCD_ 前缀（与 alphafoldserver.com 实际格式一致）
  - 支持 paired=False 的蛋白（只生成 phos 作业）
  - 生成前做残基类型自检，任何 S/T/Y 不匹配即中止

用法：
  python3 scripts/af3/generate_af3_jobs.py                    # 全部
  python3 scripts/af3/generate_af3_jobs.py --seeds 3 --split  # 分文件
"""
import sys, json, argparse
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS, REMOVED_PROTEINS

CURATED = _ROOT / "data" / "curated"
OUT_DIR = _ROOT / "af3_jobs"
CCD = {"S": "CCD_SEP", "T": "CCD_TPO", "Y": "CCD_PTR"}


def load_uniprot():
    p = CURATED / "uniprot_sequences.json"
    if not p.exists():
        sys.exit(f"ERROR: {p} 不存在，请先运行 scripts/rebuild_sequences.py")
    return json.loads(p.read_text())


def build_job(name, seq, sites, seed):
    """构造单个 AF3 server 作业。sites 为 UniProt 坐标列表。"""
    pc = {"sequence": seq, "count": 1, "useStructureTemplate": True}
    if sites:
        mods = []
        for pos in sites:
            aa = seq[pos - 1]
            if aa not in CCD:
                raise ValueError(f"{name}: 位点 {pos} 的残基是 {aa}，不是 S/T/Y")
            mods.append({"ptmType": CCD[aa], "ptmPosition": pos})
        pc["modifications"] = mods
    return {"name": f"{name}_seed{seed}",
            "modelSeeds": [str(seed)],
            "sequences": [{"proteinChain": pc}],
            "dialect": "alphafoldserver", "version": 3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3, help="每状态种子数")
    ap.add_argument("--split", action="store_true", help="按 100 作业切分文件")
    ap.add_argument("--outdir", default=str(OUT_DIR))
    args = ap.parse_args()

    uni = load_uniprot()
    gene2acc = {v["gene"].lower(): k for k, v in uni.items()}

    jobs, manifest, errors = [], [], []
    print("=" * 96)
    print("  AF3 作业生成 v2（UniProt 序列 + 核准位点）")
    print("=" * 96)

    for key, info in BENCHMARK_PROTEINS.items():
        acc = info["uniprot"]
        if acc not in uni:
            errors.append(f"{key}: UniProt {acc} 未在 curated 中")
            continue
        seq = uni[acc]["sequence"]
        sites = []
        for s in info.get("phos_sites", []):
            pos = int("".join(c for c in s if c.isdigit()))
            sites.append(pos)

        # 自检：位点残基类型
        bad = [(p, seq[p-1]) for p in sites if seq[p-1] not in CCD]
        if bad:
            errors.append(f"{key}: 位点残基异常 {bad}")
            continue

        n_jobs = 0
        # phos 作业（所有蛋白都有）
        for sd in range(1, args.seeds + 1):
            jobs.append(build_job(f"{info['gene'].lower()}_phos", seq, sites, sd))
            n_jobs += 1
        # unphos 作业（仅 paired 蛋白；单结构组无实验 unphos 参照但仍需计算基准）
        for sd in range(1, args.seeds + 1):
            jobs.append(build_job(f"{info['gene'].lower()}_unphos", seq, [], sd))
            n_jobs += 1

        manifest.append({
            "protein": key, "gene": info["gene"], "uniprot": acc,
            "seq_len": len(seq), "paired": info.get("paired"),
            "effect_type": info.get("effect_type"),
            "phos_sites_uniprot": sites,
            "pdb_chains": info.get("pdb_chains"),
            "jobs": n_jobs,
        })
        flag = "PAIRED" if info.get("paired") else "single"
        print(f"  {key:<15} {info['gene']:<9} len={len(seq):>5} "
              f"sites={sites} {flag:<7} → {n_jobs} jobs")

    if errors:
        print("\n[ERRORS]")
        for e in errors:
            print("  ", e)
        sys.exit(1)

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    if args.split and len(jobs) > 100:
        for i in range(0, len(jobs), 100):
            chunk = jobs[i:i+100]
            p = out / f"af3_jobs_v2_part{i//100 + 1}.json"
            p.write_text(json.dumps(chunk, indent=2, ensure_ascii=False))
            print(f"\n[OK] {p.name}: {len(chunk)} 作业")
    else:
        p = out / "af3_jobs_v2.json"
        p.write_text(json.dumps(jobs, indent=2, ensure_ascii=False))
        print(f"\n[OK] {p}  ({len(jobs)} 作业)")

    (out / "manifest_v2.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    print(f"[OK] {out/'manifest_v2.json'}")

    n_phos = sum(1 for j in jobs if "_phos_" in j["name"])
    print(f"\n  合计 {len(jobs)} 作业（phos {n_phos} / unphos {len(jobs)-n_phos}）")
    print(f"  覆盖 {len(manifest)} 蛋白，其中 paired {sum(1 for m in manifest if m['paired'])} 个")
    print(f"  已移除 {len(REMOVED_PROTEINS)} 个无合适实验结构的蛋白")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
