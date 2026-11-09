#!/usr/bin/env python3
"""
build_manifest.py — 冻结唯一权威基准 manifest (B1)
================================================================
Build and freeze the single authoritative benchmark manifest (blocker B1).

背景 / Background
-----------------
仓库内"纯净配对"计数存在 n=2 / n=6 / n=10 / n=11 多套互斥口径
（audit/CONFOUNDING_ANALYSIS.md、PURE_PAIR_REPORT.md、EXPANSION_REPORT.md、
audit/FIG4_HYPOTHESIS_REFUTED.md），且 data/benchmark_dataset.py 的 TOTAL_SUMMARY
字符串（"Total 10 / Paired 6"）与代码实际（18 条目 / 16 paired）严重不符。
本脚本以 data/benchmark_dataset.py 为唯一蛋白定义源，结合 audit/pair_audit.json
的 C1–C5 门禁判定与 results/af3/models/v2 的实际模型可用性，产出**唯一** manifest，
供全部下游分析与论文 Table 1 引用。

纯度判据 / Purity criteria（全部满足 = PASS）
--------------------------------------------
  P1 (=C1) 身份核验 HIT：PDB 头部含声明基因/描述关键词
  P2 (=C2) phos 态含真实磷酸化（MODRES SEP/TPO/PTR）
  P3 (=C3) phos 与 unphos 链与 UniProt 覆盖度均 ≥ 0.5
  P4 (=C4) 两态残基范围重叠 ≥ 20 且覆盖目标位点
  P5 (=C5) 两态链组成完全一致（组成纯净）

队列定义 / Cohorts
------------------
  primary (主分析队列): paired 且两态 AF3 模型齐全 → 零相关主结果
  pure    (纯度分层队列): primary ∩ P1–P5 全过
  single  : 仅有 phos 实验结构，无法做两态位移比较
  removed : 已从基准移除的条目（保留追溯）

输出 / Outputs
--------------
  results/manifest_frozen.json   机器可读权威 manifest
  audit/MANIFEST_FROZEN.md       人读报告（论文 Table 1 依据）

用法 / Usage
------------
  python3 scripts/build_manifest.py
"""
import sys, json, hashlib
from datetime import datetime, timezone
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS, REMOVED_PROTEINS  # noqa: E402

DATASET = _ROOT / "data" / "benchmark_dataset.py"
AUDIT = _ROOT / "audit" / "pair_audit.json"
CIF = _ROOT / "results" / "af3" / "models" / "v2"
OUT_JSON = _ROOT / "results" / "manifest_frozen.json"
OUT_MD = _ROOT / "audit" / "MANIFEST_FROZEN.md"

MIN_COV = 0.5          # P3 覆盖度阈值（与 audit_pairs.py 默认一致）
MIN_OVERLAP = 20       # P4 两态重叠阈值


def sha256(path):
    h = hashlib.sha256()
    h.update(Path(path).read_bytes())
    return h.hexdigest()[:16]


def af3_models(gene):
    """返回 {state: n_model_cif}；只统计 fold_*_seed*_model_*.cif（排除 template_hit）"""
    base = CIF / gene.lower()
    out = {}
    for state in ("phos", "unphos"):
        d = base / state
        n = len(list(d.glob("fold_*_seed*_model_*.cif"))) if d.exists() else 0
        out[state] = n
    return out


def main():
    audit = {r["protein"]: r for r in json.loads(AUDIT.read_text())}

    proteins = []
    for key, info in BENCHMARK_PROTEINS.items():
        gene = info["gene"]
        a = audit.get(key, {})
        pdbs = info.get("pdbs", {})
        pdb_chains = info.get("pdb_chains", {})
        sites = [int("".join(c for c in s if c.isdigit())) for s in info.get("phos_sites", [])]

        up_ids = pdbs.get("unphos", [])
        ph_ids = pdbs.get("phos", [])
        up_id = up_ids[0] if up_ids else None
        ph_id = ph_ids[0] if ph_ids else None

        # 逐 PDB 的 C1/C2/C3 明细
        def pdb_detail(pid):
            if not pid:
                return None
            e = a.get("pdbs", {}).get(pid.upper(), {})
            return {
                "pdb": pid.upper(),
                "chain": pdb_chains.get(pid, pdb_chains.get(pid.upper(), info.get("chain", "A"))),
                "title": (e.get("title") or "")[:110],
                "c1_identity": e.get("c1_identity"),
                "c2_has_real_phos": e.get("c2_has_real_phos"),
                "c2_modres_chain": e.get("c2_modres_chain", {}),
                "c3_cov_vs_uniprot": e.get("c3_cov_vs_uniprot"),
                "c3_range": e.get("c3_range"),
            }

        det_up, det_ph = pdb_detail(up_id), pdb_detail(ph_id)
        c4 = a.get("c4", {}) or {}
        overlap = c4.get("overlap_residues")
        sites_in_ph = c4.get("sites_in_phos_range", [])
        sites_in_up = c4.get("sites_in_unphos_range", [])

        af3 = af3_models(gene)
        af3_up, af3_ph = af3["unphos"], af3["phos"]
        af3_available = (af3_up > 0) and (af3_ph > 0)

        paired = bool(info.get("paired")) and bool(up_id) and bool(ph_id)

        # ---- P1–P5 ----
        p1 = bool((det_up and det_up["c1_identity"] == "HIT") and
                  (det_ph and det_ph["c1_identity"] == "HIT")) if paired else \
             bool(det_ph and det_ph["c1_identity"] == "HIT")
        p2 = bool(det_ph and det_ph["c2_has_real_phos"])
        covs = [d["c3_cov_vs_uniprot"] for d in (det_up, det_ph)
                if d and d["c3_cov_vs_uniprot"] is not None]
        p3 = bool(covs) and all(c >= MIN_COV for c in covs)
        p4 = bool(paired and overlap is not None and overlap >= MIN_OVERLAP and
                  sites and len(sites_in_ph) == len(sites))
        p5 = bool(paired and a.get("c5_pure_pair") is True)
        verdict = "PASS" if all([p1, p2, p3, p4, p5]) else a.get("verdict", "NA")

        in_primary = paired and af3_available
        is_pure = in_primary and all([p1, p2, p3, p4, p5])

        # ---- 排除原因 ----
        reason = None
        if not paired:
            reason = "single_state: 仅有 phos 实验结构"
        elif not af3_available:
            reason = f"af3_missing: AF3 模型缺失（phos={af3_ph}, unphos={af3_up}）"

        proteins.append({
            "key": key, "gene": gene, "uniprot": info["uniprot"],
            "description": info.get("description", ""),
            "effect_type": info.get("effect_type"),
            "phos_type": info.get("phos_type"),
            "sites_label": info.get("phos_sites", []),
            "sites_uniprot": sites,
            "paired": paired,
            "unphos_pdb": up_id, "phos_pdb": ph_id,
            "unphos_chain": (det_up or {}).get("chain"),
            "phos_chain": (det_ph or {}).get("chain"),
            "unphos_detail": det_up, "phos_detail": det_ph,
            "c4_overlap_residues": overlap,
            "c4_sites_in_phos": sites_in_ph, "c4_sites_in_unphos": sites_in_up,
            "c5_chains_unphos": a.get("c5_chains_unphos"),
            "c5_chains_phos": a.get("c5_chains_phos"),
            "af3_models": af3, "af3_available": af3_available,
            "criteria": {"P1_identity": p1, "P2_real_phos": p2, "P3_coverage": p3,
                         "P4_overlap": p4, "P5_composition": p5},
            "verdict": verdict,
            "in_primary_cohort": in_primary,
            "pure_pair": is_pure,
            "exclusion_reason": reason,
        })

    primary = sorted([p["key"] for p in proteins if p["in_primary_cohort"]])
    pure = sorted([p["key"] for p in proteins if p["pure_pair"]])
    single = sorted([p["key"] for p in proteins if not p["paired"]])
    fail = sorted([p["key"] for p in proteins if p["paired"] and p["verdict"] == "FAIL"])

    counts = {
        "benchmark_total": len(BENCHMARK_PROTEINS),
        "paired": sum(1 for p in proteins if p["paired"]),
        "single": len(single),
        "af3_paired_available": len(primary),
        "primary_cohort": len(primary),
        "pure_pass": sum(1 for p in proteins if p["verdict"] == "PASS"),
        "pure_cohort": len(pure),
        "warn": sum(1 for p in proteins if p["verdict"] == "WARN"),
        "fail": len(fail),
        "removed": len(REMOVED_PROTEINS),
    }

    manifest = {
        "frozen_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source": {
            "dataset": str(DATASET.relative_to(_ROOT)), "dataset_sha256": sha256(DATASET),
            "audit": str(AUDIT.relative_to(_ROOT)), "audit_sha256": sha256(AUDIT),
            "af3_models": str(CIF.relative_to(_ROOT)),
        },
        "criteria": {"min_coverage": MIN_COV, "min_overlap": MIN_OVERLAP,
                     "P1": "identity HIT", "P2": "phos has SEP/TPO/PTR",
                     "P3": "coverage >= 0.5 both states",
                     "P4": "overlap >= 20 residues and covers target site",
                     "P5": "identical chain composition"},
        "cohorts": {"primary": primary, "pure": pure, "single": single,
                    "fail": fail, "removed": sorted(REMOVED_PROTEINS.keys())},
        "counts": counts,
        "proteins": proteins,
    }

    OUT_JSON.write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    _write_md(manifest)

    print("=" * 92)
    print("  基准 manifest 冻结 (B1)")
    print("=" * 92)
    print(f"  基准总数 {counts['benchmark_total']} | 成对 {counts['paired']} | 单结构 {counts['single']} "
          f"| 已移除 {counts['removed']}")
    print(f"  PASS {counts['pure_pass']} | WARN {counts['warn']} | FAIL {counts['fail']}")
    print(f"  主分析队列 (paired + AF3 齐全) n={counts['primary_cohort']}")
    print(f"  纯净队列 (P1–P5 全过)          n={counts['pure_cohort']}")
    print("-" * 92)
    print(f"  {'蛋白':<24}{'P1':>3}{'P2':>3}{'P3':>3}{'P4':>3}{'P5':>3}  {'判定':<6}{'主队列':<7}{'纯净'}")
    for p in proteins:
        c = p["criteria"]
        mark = lambda b: "✓" if b else "·"
        print(f"  {p['key']:<24}{mark(c['P1_identity']):>3}{mark(c['P2_real_phos']):>3}"
              f"{mark(c['P3_coverage']):>3}{mark(c['P4_overlap']):>3}{mark(c['P5_composition']):>3}  "
              f"{p['verdict']:<6}{'✓' if p['in_primary_cohort'] else '·':<7}"
              f"{'✓' if p['pure_pair'] else '·'}")
    print(f"\n[OK] {OUT_JSON}")
    print(f"[OK] {OUT_MD}")
    return manifest


def _write_md(m):
    c = m["counts"]
    lines = [
        "# 基准 Manifest（冻结版）", "",
        f"冻结时间：{m['frozen_at']}　来源：`{m['source']['dataset']}` "
        f"(sha256 `{m['source']['dataset_sha256']}`)、`{m['source']['audit']}`", "",
        "> **本文件为基准唯一权威记录。** 论文 Table 1 与全部下游分析必须引用此 manifest；",
        "> 历史文档（CONFOUNDING_ANALYSIS / PURE_PAIR_REPORT / EXPANSION_REPORT）中的 n=2/6/10/11 ",
        "> 口径均以本表为准作废。", "",
        "## 计数 / Counts", "",
        "| 项 | 数量 |", "|---|---|",
        f"| 基准总数 benchmark_total | {c['benchmark_total']} |",
        f"| 成对 paired | {c['paired']} |",
        f"| 单结构 single | {c['single']} |",
        f"| 主分析队列 primary（paired + AF3 齐全） | {c['primary_cohort']} |",
        f"| 纯净队列 pure（P1–P5 全过） | {c['pure_cohort']} |",
        f"| 判定 PASS / WARN / FAIL | {c['pure_pass']} / {c['warn']} / {c['fail']} |",
        f"| 已移除 removed | {c['removed']} |", "",
        "## 判据 / Criteria", "",
        f"- P3 覆盖度阈值 = {m['criteria']['min_coverage']}",
        f"- P4 两态重叠阈值 = {m['criteria']['min_overlap']} 残基",
        "- P1 身份核验；P2 phos 态含 SEP/TPO/PTR；P5 两态链组成完全一致", "",
        "## 逐蛋白 / Per-protein", "",
        "| 蛋白 | 基因 | UniProt | 位点 | phos PDB | unphos PDB | 判定 | P1 | P2 | P3 | P4 | P5 | 主队列 | 纯净 |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for p in m["proteins"]:
        cr = p["criteria"]
        mk = lambda b: "✓" if b else "·"
        lines.append(
            f"| {p['key']} | {p['gene']} | {p['uniprot']} | {','.join(map(str, p['sites_uniprot']))} "
            f"| {p['phos_pdb'] or '—'} | {p['unphos_pdb'] or '—'} | {p['verdict']} "
            f"| {mk(cr['P1_identity'])} | {mk(cr['P2_real_phos'])} | {mk(cr['P3_coverage'])} "
            f"| {mk(cr['P4_overlap'])} | {mk(cr['P5_composition'])} "
            f"| {'✓' if p['in_primary_cohort'] else '·'} | {'✓' if p['pure_pair'] else '·'} |")
    lines += ["", "## 队列成员 / Cohort membership", "",
              f"- **primary (n={c['primary_cohort']})**：{', '.join(m['cohorts']['primary'])}",
              f"- **pure (n={c['pure_cohort']})**：{', '.join(m['cohorts']['pure'])}",
              f"- **single (n={c['single']})**：{', '.join(m['cohorts']['single'])}",
              f"- **fail (n={c['fail']})**：{', '.join(m['cohorts']['fail'])}",
              f"- **removed (n={c['removed']})**：{', '.join(m['cohorts']['removed'])}", ""]
    OUT_MD.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
