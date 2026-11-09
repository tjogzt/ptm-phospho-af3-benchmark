#!/usr/bin/env python3
"""
make_tables.py — 由冻结产物生成论文 Table 1 / Table 2 / Table S1
==================================================================
Build the paper's Table 1 / Table 2 / Table S1 from frozen artifacts.

中文（Chinese）
---------------
本脚本**只读取冻结产物**，不做任何重新计算，保证表内数字与正文/图完全同源：
  - results/manifest_frozen.json          （Table 1：基准与纯度门禁 P1–P5）
  - results/analysis/af3_vs_pdb_summary.json（Table 2：逐蛋白指标；Table S1：队列统计）
  - results/analysis/esmfold_vs_pdb_summary.json（存在时才并入对比线）

输出（同时给出 CSV 与 booktabs LaTeX 片段）：
  results/tables/Table1_benchmark.{csv,tex}
  results/tables/Table2_metrics.{csv,tex}
  results/tables/TableS1_stats.{csv,tex}

English
-------
Reads ONLY frozen artifacts (no recomputation) so the tables are numerically
identical to the figures and the manuscript body. Emits both CSV and
booktabs LaTeX fragments.

用法 / Usage：
    python3 scripts/analysis/make_tables.py
"""
import json
import csv
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent.parent
ANA = _ROOT / "results" / "analysis"
TAB = _ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

MANIFEST = _ROOT / "results" / "manifest_frozen.json"

# 效应类型 → 简洁缩写（表内英文标注）/ Effect-type abbreviations (English labels)
EFFECT_SHORT = {
    "local_rearrangement": "local",
    "allosteric": "allost.",
    "disorder_order": "disorder",
}
# 磷酸化残基类型 → 单字母 + p 前缀 / Phospho-residue type to one-letter with p prefix
PHOS_MAP = {"SEP": "pS", "TPO": "pT", "PTR": "pY", "SEP_TPO": "pS/pT",
            "TPO_PTR": "pT/pY", "SEP_PTR": "pS/pY"}


def _phos_label(ptype):
    """把 manifest 的 phos_type（如 TPO_PTR）转成 pT/pY 形式。"""
    if ptype in PHOS_MAP:
        return PHOS_MAP[ptype]
    return "/".join(PHOS_MAP.get(p, p) for p in str(ptype).split("_"))


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"[OK] {path}")


def write_tex(path, caption, label, header, rows, colspec, resize=False):
    """生成可直接 \\input 的 booktabs 表。resize=True 时按 \\textwidth 等比缩放以适配页宽。"""
    lines = [
        "\\begin{table}[htbp]",
        "\\centering",
    ]
    if resize:
        lines.append("\\resizebox{\\textwidth}{!}{%")
    lines += [
        f"\\begin{{tabular}}{{{colspec}}}",
        "\\toprule",
        " & ".join(header) + " \\\\",
        "\\midrule",
    ]
    for r in rows:
        lines.append(" & ".join(str(c) for c in r) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    if resize:
        lines.append("}")
    lines += [
        f"\\caption{{{caption}}}",
        f"\\label{{{label}}}",
        "\\end{table}",
        "",
    ]
    Path(path).write_text("\n".join(lines))
    print(f"[OK] {path}")


def main():
    mf = json.loads(MANIFEST.read_text())
    af3 = json.loads((ANA / "af3_vs_pdb_summary.json").read_text())
    # 注意：af3 summary 的 "protein" 键 = manifest 的 key（如 ERK2 / c-Src_kinase），
    # 而非 gene（如 MAPK1 / SRC）；须按 key 查表，否则会静默丢蛋白。
    # NOTE: the af3 summary keys on the manifest `key` (e.g. ERK2), not `gene` (MAPK1).
    af3_by = {r["protein"]: r for r in af3["per_protein"]}

    # 需要出现在 Table 1/2 的队列 = primary（15）+ 被排除者（供审计透明）
    primary = list(mf["cohorts"]["primary"])
    excluded = [p["key"] for p in mf["proteins"]
                if not p["in_primary_cohort"]]

    # ---------------- Table 1：基准与纯度门禁 ----------------
    t1_h = ["Protein", "UniProt", "Effect", "PTM", "Site(s)",
            "Unphos PDB", "Phos PDB", "P1", "P2", "P3", "P4", "P5",
            "Verdict", "Pure", "In cohort"]
    t1_rows = []
    for p in mf["proteins"]:
        c = p["criteria"]
        yn = lambda b: "Y" if b else ("N" if b is False else "--")
        sites = "/".join(p.get("sites_label") or []) if p.get("sites_label") else "--"
        # 位点后附磷酸化类型，如 "Thr160 (pT)"
        sites_txt = f"{sites} ({_phos_label(p['phos_type'])})" if sites != "--" else "--"
        t1_rows.append([
            p["gene"], p["uniprot"],
            EFFECT_SHORT.get(p["effect_type"], p["effect_type"]),
            _phos_label(p["phos_type"]),
            sites_txt,
            p.get("unphos_pdb") or "--",
            p.get("phos_pdb") or "--",
            yn(c.get("P1_identity")), yn(c.get("P2_real_phos")),
            yn(c.get("P3_coverage")), yn(c.get("P4_overlap")),
            yn(c.get("P5_composition")),
            p["verdict"], "Y" if p["pure_pair"] else "N",
            "Y" if p["in_primary_cohort"] else "N",
        ])
    write_csv(TAB / "Table1_benchmark.csv", t1_h, t1_rows)
    write_tex(TAB / "Table1_benchmark.tex",
              "Benchmark cohort and pairing-purity gate (P1--P5). "
              "P1 identity, P2 real phospho-residue (MODRES), "
              "P3 coverage $\\geq$0.5 in both states, P4 $\\geq$20-residue "
              "overlap covering the target site, P5 identical chain composition.",
              "tab:benchmark", t1_h, t1_rows,
              "lllll lcc ccccc cc", resize=True)

    # ---------------- Table 2：逐蛋白指标 ----------------
    has_esm = (ANA / "esmfold_vs_pdb_summary.json").exists()
    esm_by = {}
    if has_esm:
        esm = json.loads((ANA / "esmfold_vs_pdb_summary.json").read_text())
        esm_by = {r["protein"]: r for r in esm["per_protein"]}

    t2_h = ["Protein", "UniProt", "$n_{common}$",
            "exp disp (\\AA)", "AF3 disp (\\AA)", "AF3 ratio",
            "AF3 cos", "AF3 $r$"]
    if has_esm:
        t2_h += ["ESM disp (\\AA)", "ESM ratio", "ESM cos", "ESM $r$"]
    t2_h += ["Pure"]

    t2_rows = []
    for key in primary:
        p = next(x for x in mf["proteins"] if x["key"] == key)
        a = af3_by.get(key)
        if a is None:
            continue
        exp_d = a["mean_disp_exp"]
        pred_d = a["mean_disp_pred"]
        ratio = (pred_d / exp_d) if exp_d else float("nan")
        cos = a.get("mean_cosine")
        cos_txt = f"{cos:+.2f}" if cos is not None else "--"
        row = [
            p["gene"], p["uniprot"], a.get("n_common", "--"),
            f"{exp_d:.2f}", f"{pred_d:.2f}", f"{ratio:.2f}",
            cos_txt, f"{a.get('pearson_magnitude', float('nan')):+.2f}",
        ]
        if has_esm:
            e = esm_by.get(key)
            if e:
                ed = e["mean_disp_pred"]
                er = (ed / exp_d) if exp_d else float("nan")
                ecos = e.get("mean_cosine")
                row += [f"{ed:.2f}", f"{er:.2f}",
                        (f"{ecos:+.2f}" if ecos is not None else "--"),
                        f"{e.get('pearson_magnitude', float('nan')):+.2f}"]
            else:
                row += ["--", "--", "--", "--"]
        row += ["Y" if p["pure_pair"] else "N"]
        t2_rows.append(row)

    write_csv(TAB / "Table2_metrics.csv", t2_h, t2_rows)
    write_tex(TAB / "Table2_metrics.tex",
              "Per-protein comparison of experimental versus predicted "
              "phospho-induced C$\\alpha$ displacement. Displacements are "
              "mean per-residue magnitudes (\\AA). Ratio $=$ predicted/experimental. "
              "cos $=$ mean direction cosine over residues moving $>$1\\,\\AA. "
              "$r$ $=$ Pearson $r$ of the per-residue magnitude profiles.",
              "tab:metrics", t2_h, t2_rows,
              ("ll c cccccc ccc c" if has_esm else "ll c ccccc c"), resize=True)

    # ---------------- Table S1：队列统计 ----------------
    cs = af3["cohort_stats"]
    s1_rows = [
        ["$n$ proteins", cs["n_proteins"]],
        ["Pearson $r$ (pred vs exp)", f"{cs['pearson_r']:+.4f}"],
        ["Pearson 95\\% CI", f"[{cs['pearson_CI95'][0]:+.4f}, {cs['pearson_CI95'][1]:+.4f}]"],
        ["slope (pred per exp, \\AA/\\AA)", f"{cs['slope_pred_per_exp']:+.4f}"],
        ["Spearman $\\rho$", f"{cs['spearman_rho']:+.4f}"],
        ["Spearman $p$", f"{cs['spearman_p']:.4f}"],
        ["mean predicted disp (\\AA)", f"{cs['mean_pred']:.3f}"],
        ["SD predicted disp (\\AA)", f"{cs['sd_pred']:.3f}"],
        ["experimental disp range (\\AA)",
         f"[{cs['range_exp'][0]:.2f}, {cs['range_exp'][1]:.2f}]"],
        ["max leverage (removed)", cs["max_leverage"]["removed"]],
        ["$r$ without max-leverage point", f"{cs['max_leverage']['r_without']:+.4f}"],
    ]
    if has_esm:
        ecs = esm["cohort_stats"]
        s1_rows += [
            ["ESMFold Pearson $r$", f"{ecs['pearson_r']:+.4f}"],
            ["ESMFold Pearson 95\\% CI",
             f"[{ecs['pearson_CI95'][0]:+.4f}, {ecs['pearson_CI95'][1]:+.4f}]"],
            ["ESMFold mean predicted disp (\\AA)", f"{ecs['mean_pred']:.3f}"],
            ["ESMFold max leverage (removed)", ecs["max_leverage"]["removed"]],
            ["ESMFold $r$ without max-leverage point",
             f"{ecs['max_leverage']['r_without']:+.4f}"],
        ]
    write_csv(TAB / "TableS1_stats.csv",
              ["Quantity", "Value"],
              [[a, b] for a, b in s1_rows])
    write_tex(TAB / "TableS1_stats.tex",
              "Cohort-level statistics for the AF3 (and, where available, "
              "ESMFold) predicted-versus-experimental displacement comparison. "
              "Bootstrap CIs resample proteins (the independent unit), not models.",
              "tab:stats", ["Quantity", "Value"], s1_rows, "ll")

    print(f"\n[done] has_esm={has_esm}; tables in {TAB}")


if __name__ == "__main__":
    main()
