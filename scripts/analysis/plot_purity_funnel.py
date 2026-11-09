#!/usr/bin/env python3
"""
plot_purity_funnel.py — Fig 1（Hero 图）：配对纯度漏斗 + 纯度分层对比
====================================================================
Fig 1 (hero) — pairing-purity funnel and purity-stratified comparison.

中文（Chinese）
---------------
左panel：由候选池经 C1–C5（P1–P5）门禁收敛到可用纯配对的漏斗，数字全部取自
         冻结 manifest（counts / cohorts），不重新推导。
右panel：把同一份**冻结**逐蛋白结果按配对纯度分层（pure vs impure vs all），
         比较三类表观指标（幅度比中位数、方向余弦中位数、幅度 Pearson r）。
         用于展示「配对纯度如何改变结论」。

English
-------
Left panel: the candidate pool funnelled through the P1--P5 purity gates down to
the usable pure pairs; all numbers are read verbatim from the frozen manifest.
Right panel: the same frozen per-protein results stratified by pairing purity
(pure vs impure vs all) across three apparent metrics (median magnitude ratio,
median direction cosine, Pearson r of the magnitude profiles), illustrating how
pairing purity changes the conclusion.

输入（全部为冻结产物）/ Inputs (frozen artifacts only):
  results/manifest_frozen.json
  results/analysis/af3_vs_pdb_summary.json
输出 / Output:
  results/figures/Fig1_purity_funnel.pdf / .png
用法 / Usage:
  python3 scripts/analysis/plot_purity_funnel.py
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ROOT = Path(__file__).resolve().parent.parent.parent
ANA = _ROOT / "results" / "analysis"
FIG = _ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
                     "axes.labelsize": 7.5, "axes.titlesize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42})

mf = json.loads((_ROOT / "results" / "manifest_frozen.json").read_text())
cnt = mf["counts"]
PURE = set(mf["cohorts"]["pure"])

af3 = json.loads((ANA / "af3_vs_pdb_summary.json").read_text())
rows = [r for r in af3["per_protein"] if not r.get("same_structure")]


def _sub(members):
    return [r for r in rows if r["protein"] in members]


def _stats(sub):
    """纯描述统计（中位数 + Pearson r）。"""
    if not sub:
        return None
    exp = np.array([r["mean_disp_exp"] for r in sub], float)
    prd = np.array([r["mean_disp_pred"] for r in sub], float)
    ratio = float(np.median(prd / exp))
    cs = [r["mean_cosine"] for r in sub if r.get("mean_cosine") is not None]
    cos = float(np.median(cs)) if cs else float("nan")
    r = float(np.corrcoef(exp, prd)[0, 1]) if len(sub) > 2 else float("nan")
    return {"n": len(sub), "ratio": ratio, "cos": cos, "r": r}


pure_rows = _sub(PURE)
impure_rows = _sub({r["protein"] for r in rows} - PURE)
S_P, S_I, S_A = _stats(pure_rows), _stats(impure_rows), _stats(rows)

C_PURE, C_IMP, C_ALL = "#2ca02c", "#d62728", "#7f7f7f"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0),
                               gridspec_kw={"width_ratios": [1.15, 1.0]})

# ================= Panel a：纯度漏斗 =================
# 阶段与计数全部来自冻结 manifest
stages = [
    ("Curated candidate pool",
     cnt["benchmark_total"] + cnt["removed"]),
    ("Two-state pairs (unphos + phos)", cnt["paired"]),
    ("AF3 models for both states", cnt["primary_cohort"]),
    ("Passed purity gate P1\u2013P5", cnt["pure_cohort"]),
]
n0 = stages[0][1]
ys = np.arange(len(stages))[::-1]
for i, (lab, val) in enumerate(stages):
    y = ys[i]
    ax1.barh(y, val, height=0.62,
             color=plt.cm.Blues(0.30 + 0.55 * (val / n0)),
             edgecolor="black", linewidth=0.5)
    ax1.text(val + n0 * 0.012, y, f"{val}", va="center", ha="left", fontsize=7.2)
    ax1.text(0, y + 0.46, lab, va="bottom", ha="left", fontsize=6.6)

ax1.set_yticks([])
ax1.set_xlim(0, n0 * 1.14)
ax1.set_xlabel("Number of benchmark entries")
ax1.set_title("a   Pairing-purity funnel", loc="left", fontweight="bold")
ax1.spines[["top", "right", "left"]].set_visible(False)
ax1.tick_params(left=False)
# 门禁判据注释
ax1.text(n0 * 0.02, -1.05,
         "P1 identity  \u00b7  P2 MODRES phospho-residue  \u00b7  P3 coverage "
         "\u22650.5  \u00b7  P4 \u226520-res overlap  \u00b7  P5 identical chains",
         fontsize=5.9, style="italic", color="0.35")

# ================= Panel b：纯度分层对比 =================
metrics = [("Magnitude\nratio", "ratio"), ("Direction\ncosine", "cos"),
           ("Magnitude\n$r$", "r")]
groups = [("Pure", S_P, C_PURE), ("Impure", S_I, C_IMP), ("All", S_A, C_ALL)]
nb = len(metrics)
w = 0.26
x0 = np.arange(nb)
for gi, (gname, gs, col) in enumerate(groups):
    xs = x0 + (gi - 1) * w
    vals = [gs[m[1]] for m in metrics]
    ax2.bar(xs, vals, width=w, color=col, edgecolor="black", linewidth=0.5,
            label=f"{gname} (n={gs['n']})")
ax2.axhline(0, color="0.5", lw=0.6)
ax2.set_xticks(x0)
ax2.set_xticklabels([m[0] for m in metrics])
ax2.set_ylabel("Value")
ax2.set_title("b   Stratified by pairing purity", loc="left", fontweight="bold")
ax2.legend(frameon=False, loc="upper right", handletextpad=0.4,
           borderpad=0.2, labelspacing=0.3)
ax2.spines[["top", "right"]].set_visible(False)

fig.tight_layout(pad=0.9, w_pad=1.8)
fig.savefig(FIG / "Fig1_purity_funnel.pdf", bbox_inches="tight")
fig.savefig(FIG / "Fig1_purity_funnel.png", bbox_inches="tight", dpi=300)
plt.close(fig)

print("[OK]", FIG / "Fig1_purity_funnel.pdf")
print(f"  pure   (n={S_P['n']}) ratio={S_P['ratio']:.3f} cos={S_P['cos']:.3f} r={S_P['r']:+.3f}")
print(f"  impure (n={S_I['n']}) ratio={S_I['ratio']:.3f} cos={S_I['cos']:.3f} r={S_I['r']:+.3f}")
print(f"  all    (n={S_A['n']}) ratio={S_A['ratio']:.3f} cos={S_A['cos']:.3f} r={S_A['r']:+.3f}")
