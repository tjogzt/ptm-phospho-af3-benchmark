#!/usr/bin/env python3
"""
analyze_af3_ptm_effect.py — P0-4 第一版分析：AF3 磷酸化结构效应量化
======================================================================
输入：results/af3/af3_summary.csv（由 download_af3_results.py 幂等生成）
输出：
  results/figures/Fig1_af3_ptm_effect.pdf / .png   双面板主图
  results/figures/Fig1_af3_ptm_effect_stats.csv    逐蛋白统计表
  results/analysis_summary.json                     汇总指标
规范：英文标注、最小字号 7pt、PDF 矢量输出、配色全局一致、目标尺寸直出。

用法：python scripts/analysis/analyze_af3_ptm_effect.py
"""
import json, csv
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_ROOT = Path(__file__).resolve().parent.parent.parent
AF3 = _ROOT / "results" / "af3"
FIG = _ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# SCI 图表规范：统一配色与字号
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7,
    "axes.labelsize": 7.5, "axes.titlesize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
})
COLORS = {
    "local_rearrangement": "#1f77b4",
    "allosteric":          "#d62728",
    "disorder_order":      "#2ca02c",
    "negative_control":    "#7f7f7f",
}
LABEL = {
    "local_rearrangement": "Local rearrangement",
    "allosteric":          "Allosteric",
    "disorder_order":      "Disorder$\\rightarrow$order",
    "negative_control":    "Negative control",
}
# 短名（图内标注）
SHORT = {
    "CDK2": "CDK2", "ERK2": "ERK2", "PKA_Calpha": "PKA", "Aurora_A": "AURKA",
    "CDK7": "CDK7", "P38alpha": "p38", "IRF3": "IRF3", "CFTR_R_domain": "CFTR",
    "Hsp27": "HSP27", "AMPK_alpha2": "AMPK", "PhK_gamma": "PhK",
    "c-Src_kinase": "Src", "Glycogen_Phosphorylase": "PYGM", "Ubiquitin": "UBB",
}


def welch(a, b):
    """Welch t 检验（返回 t, 近似 p；无 scipy 依赖）。"""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    va, vb = a.var(ddof=1)/len(a), b.var(ddof=1)/len(b)
    se = np.sqrt(va+vb)
    if se == 0: return np.nan, np.nan
    t = (a.mean()-b.mean())/se
    df = (va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1))
    # t 分布生存函数（正态近似，n<10 时偏保守提示）
    from math import erf, sqrt
    p = 2*(1-0.5*(1+erf(abs(t)/sqrt(2))))
    return float(t), float(p)


def main():
    rows = list(csv.DictReader(open(AF3/"af3_summary.csv")))
    agg = defaultdict(lambda: {"phos": {"plddt": [], "ptm": []},
                               "unphos": {"plddt": [], "ptm": []},
                               "effect": None, "seqlen": None})
    for r in rows:
        p = r["protein"]
        if p.startswith("REMOVED:") or p == "UNKNOWN":
            continue                      # 已移除/未识别条目不参与统计
        if r.get("data_version") == "v1":
            continue                      # 仅使用位点经核准的 v2 数据
        st = r["state"]
        agg[p]["effect"] = r["effect_type"]
        agg[p]["seqlen"] = int(r["seq_len"]) if r["seq_len"] else None
        if r["plddt_mean"]:
            agg[p][st]["plddt"].append(float(r["plddt_mean"]))
        if r["ptm_mean"]:
            agg[p][st]["ptm"].append(float(r["ptm_mean"]))

    stats = []
    for p, d in sorted(agg.items()):
        if not d["phos"]["plddt"] or not d["unphos"]["plddt"]:
            continue
        dp = np.mean(d["phos"]["plddt"]) - np.mean(d["unphos"]["plddt"])
        t, pv = welch(d["phos"]["plddt"], d["unphos"]["plddt"])
        dtm = (np.mean(d["phos"]["ptm"]) - np.mean(d["unphos"]["ptm"])
               if d["phos"]["ptm"] and d["unphos"]["ptm"] else np.nan)
        stats.append({
            "protein": p, "short": SHORT.get(p, p), "effect": d["effect"],
            "seq_len": d["seqlen"],
            "n_phos": len(d["phos"]["plddt"]), "n_unphos": len(d["unphos"]["plddt"]),
            "plddt_phos": round(float(np.mean(d["phos"]["plddt"])), 2),
            "plddt_unphos": round(float(np.mean(d["unphos"]["plddt"])), 2),
            "delta_plddt": round(float(dp), 2),
            "ptm_phos": round(float(np.mean(d["phos"]["ptm"])), 3) if d["phos"]["ptm"] else None,
            "ptm_unphos": round(float(np.mean(d["unphos"]["ptm"])), 3) if d["unphos"]["ptm"] else None,
            "delta_ptm": round(float(dtm), 3) if not np.isnan(dtm) else None,
            "t": None if np.isnan(t) else round(t, 2),
            "p_approx": None if np.isnan(pv) else round(pv, 4),
        })

    with open(FIG/"Fig1_af3_ptm_effect_stats.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(stats[0].keys()))
        w.writeheader(); w.writerows(stats)

    # ================= 绘图 =================
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.1, 2.9))

    # --- Panel A: ΔpLDDT vs ΔpTM 散点 ---
    # 手工微调标签偏移，避免重叠
    LABEL_OFF = {"IRF3": (5, -2), "UBB": (-16, 2), "p38": (4, -6),
                 "HSP27": (-22, 1), "CFTR": (4, -5), "PhK": (4, -6),
                 "CDK7": (4, -6), "PYGM": (4, 2), "AMPK": (4, 2)}
    for e in ["local_rearrangement", "allosteric", "disorder_order", "negative_control"]:
        sub = [s for s in stats if s["effect"] == e and s["delta_ptm"] is not None]
        if not sub: continue
        axA.scatter([s["delta_ptm"] for s in sub], [s["delta_plddt"] for s in sub],
                    s=26, color=COLORS[e], edgecolor="black", linewidth=0.4,
                    label=LABEL[e], zorder=3)
        for s in sub:
            dx, dy = LABEL_OFF.get(s["short"], (3, 3))
            axA.annotate(s["short"], (s["delta_ptm"], s["delta_plddt"]),
                         textcoords="offset points", xytext=(dx, dy), fontsize=5.6,
                         color=COLORS[e], zorder=4)
    axA.axhline(0, color="0.6", lw=0.5, ls="--", zorder=1)
    axA.axvline(0, color="0.6", lw=0.5, ls="--", zorder=1)
    axA.set_xlabel(r"$\Delta$pTM  (phos $-$ unphos)")
    axA.set_ylabel(r"$\Delta$pLDDT  (phos $-$ unphos)")
    axA.set_title("a   Confidence shift upon PTM introduction", loc="left", fontweight="bold")
    axA.legend(frameon=False, loc="upper left", handletextpad=0.4,
           borderpad=0.2, labelspacing=0.3)
    axA.spines[["top", "right"]].set_visible(False)

    # --- Panel B: 按效应类型分组的 ΔpLDDT ---
    # 仅保留至少有 1 个蛋白的效应类型（v2 数据只覆盖成对结构组）
    order = [e for e in ["local_rearrangement", "allosteric", "disorder_order", "negative_control"]
             if any(s["effect"] == e for s in stats)]
    xs, ys, cols, labels = [], [], [], []
    for i, e in enumerate(order):
        sub = sorted([s for s in stats if s["effect"] == e],
                     key=lambda s: s["delta_plddt"])
        for j, s in enumerate(sub):
            xs.append(i + (j-(len(sub)-1)/2)*0.16)
            ys.append(s["delta_plddt"]); cols.append(COLORS[e])
            labels.append(s["short"])
    axB.scatter(xs, ys, s=30, c=cols, edgecolor="black", linewidth=0.4, zorder=3)
    LAB_OFF_B = {"CFTR": (5, -1), "HSP27": (-24, 2), "IRF3": (5, -3),
                 "UBB": (5, -3), "p38": (5, -3)}
    for x, y, lb in zip(xs, ys, labels):
        dx, dy = LAB_OFF_B.get(lb, (4, 2))
        axB.annotate(lb, (x, y), textcoords="offset points", xytext=(dx, dy), fontsize=5.2)
    # 组均值横线
    for i, e in enumerate(order):
        sub = [s["delta_plddt"] for s in stats if s["effect"] == e]
        if sub:
            axB.plot([i-0.34, i+0.34], [np.mean(sub)]*2, color=COLORS[e], lw=1.6, zorder=2)
    axB.axhline(0, color="0.6", lw=0.5, ls="--", zorder=1)
    axB.set_xticks(range(len(order)))
    axB.set_xticklabels([LABEL[e].replace("Disorder$\\rightarrow$order", "Disorder→order")
                         for e in order], rotation=18, ha="right")
    axB.set_ylabel(r"$\Delta$pLDDT per protein")
    axB.set_title("b   Grouped by structural effect type", loc="left", fontweight="bold")
    axB.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(pad=0.8, w_pad=1.6)
    fig.savefig(FIG/"Fig1_af3_ptm_effect.pdf", bbox_inches="tight")
    fig.savefig(FIG/"Fig1_af3_ptm_effect.png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    # ================= 汇总 =================
    summary = {
        "n_proteins": len(stats),
        "n_jobs": len(rows),
        "by_effect": {},
        "correlation_delta_plddt_vs_ptm": None,
        "caveats": [
            "seeds 为 AF3 采样随机种子，非独立生物学重复；p 值仅反映采样噪声",
            "UBB 构建体为 685 aa 多聚泛素前体（非 76 aa 成熟泛素），负对照表述待修正",
            "CFTR 为全长 1480 aa（非 R 域片段），disorder_order 标签需重述",
            "IRF3 phospho 参考结构为 S→D  phosphomimetic（非真实磷酸化）",
        ],
    }
    for e in order:
        sub = [s["delta_plddt"] for s in stats if s["effect"] == e]
        if sub:
            summary["by_effect"][e] = {
                "n": len(sub), "mean_delta_plddt": round(float(np.mean(sub)), 2),
                "sd": round(float(np.std(sub, ddof=1)), 2) if len(sub) > 1 else None,
                "proteins": [s["short"] for s in stats if s["effect"] == e],
            }
    valid = [(s["delta_ptm"], s["delta_plddt"]) for s in stats
             if s["delta_ptm"] is not None]
    if len(valid) > 2:
        a = np.array(valid)
        summary["correlation_delta_plddt_vs_ptm"] = {
            "pearson_r": round(float(np.corrcoef(a[:, 0], a[:, 1])[0, 1]), 3),
            "n": len(valid),
        }
    (_ROOT/"results"/"analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False))

    print(f"{'='*74}")
    print(f"  分析完成：{summary['n_proteins']} 蛋白 / {summary['n_jobs']} 作业")
    print(f"{'='*74}")
    for e in order:
        if e in summary["by_effect"]:
            d = summary["by_effect"][e]
            print(f"  {LABEL[e]:28s} n={d['n']}  ΔpLDDT={d['mean_delta_plddt']:+.2f}"
                  f" ± {d['sd'] if d['sd'] is not None else '-'}")
    if summary["correlation_delta_plddt_vs_ptm"]:
        c = summary["correlation_delta_plddt_vs_ptm"]
        print(f"\n  ΔpLDDT vs ΔpTM 相关：r = {c['pearson_r']:+.3f} (n={c['n']})")
    print(f"\n  图：{FIG}/Fig1_af3_ptm_effect.pdf")
    print(f"  表：{FIG}/Fig1_af3_ptm_effect_stats.csv")
    return 0


if __name__ == "__main__":
    main()
