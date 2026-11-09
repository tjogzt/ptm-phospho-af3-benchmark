#!/usr/bin/env python3
"""
plot_af3_vs_pdb.py — Fig2：预测位移 vs PDB 实验位移（AF3 主图，含 ESMFold 双线）
======================================================================
Fig 2 是论文主结果图（C3/C4 + B7 ESMFold 对比线）。

Panel a：实验位移 vs 预测位移散点（1:1 参考线）。AF3 与 ESMFold 两套预测各一组
         点；点色按配对纯度（manifest 的 pure 队列）区分。
Panel b：方向余弦与幅度相关 r 的逐蛋白条形对比（两套预测并排）。

数据来源（均来自冻结产物，禁止硬编码）：
  results/analysis/af3_vs_pdb_summary.json
  results/analysis/esmfold_vs_pdb_summary.json   （存在时才叠加）
  results/manifest_frozen.json                    （pure 队列定义）

输出：results/figures/Fig2_af3_vs_pdb.pdf / .png
用法：python3 scripts/analysis/plot_af3_vs_pdb.py
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ROOT = Path(__file__).resolve().parent.parent.parent
ANA = _ROOT / "results" / "analysis"
FIG = _ROOT / "results" / "figures"; FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 7.5,
                     "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 6.5, "axes.linewidth": 0.6, "pdf.fonttype": 42})

# pure 队列从冻结 manifest 读取（不再硬编码）
_mf = json.loads((_ROOT / "results" / "manifest_frozen.json").read_text())
PURE_KEYS = set(_mf["cohorts"]["pure"])

af3 = json.loads((ANA / "af3_vs_pdb_summary.json").read_text())
rows = [r for r in af3["per_protein"]
        if not r.get("same_structure") and r.get("mean_disp_exp") is not None]
for r in rows:
    r["is_pure"] = r.get("protein") in PURE_KEYS or r.get("pure_pair")

esm = None
ESM_P = ANA / "esmfold_vs_pdb_summary.json"
if ESM_P.exists():
    e = json.loads(ESM_P.read_text())
    esm = [r for r in e["per_protein"]
           if not r.get("same_structure") and r.get("mean_disp_exp") is not None]
    for r in esm:
        r["is_pure"] = r.get("protein") in PURE_KEYS or r.get("pure_pair")

C_AF3, C_ESM = "#1f77b4", "#d62728"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2))

# ---------------- Panel a：幅度散点（AF3 + 可选 ESMFold） ----------------
xs = [r["mean_disp_exp"] for r in rows]
ys = [r["mean_disp_pred"] for r in rows]
allv = xs + ys
if esm:
    allv += [r["mean_disp_exp"] for r in esm] + [r["mean_disp_pred"] for r in esm]
lim = max(allv) * 1.28
ax1.plot([0, lim], [0, lim], ls="--", lw=0.7, color="0.55", zorder=1, label="1:1")

for r in rows:
    ax1.scatter(r["mean_disp_exp"], r["mean_disp_pred"], s=44, c=C_AF3,
                edgecolor="black", linewidth=0.5, zorder=3)
if esm:
    for r in esm:
        ax1.scatter(r["mean_disp_exp"], r["mean_disp_pred"], s=40, c=C_ESM,
                    marker="^", edgecolor="black", linewidth=0.5, zorder=3)

ax1.scatter([], [], s=44, c=C_AF3, edgecolor="black", linewidth=0.5, label="AF3")
if esm:
    ax1.scatter([], [], s=40, c=C_ESM, marker="^", edgecolor="black",
                linewidth=0.5, label="ESMFold (phosphomimetic)")
ax1.scatter([], [], s=44, facecolor="none", edgecolor="0.35", linewidth=0.8,
            label="Open marker = pure pair")

# 纯配对用空心描边强调
for r in rows:
    if r.get("is_pure"):
        ax1.scatter(r["mean_disp_exp"], r["mean_disp_pred"], s=44, facecolor="none",
                    edgecolor="0.35", linewidth=0.7, zorder=4)

OFF = {"CDK7": (-20, -11), "AURKA": (-8, 6), "CDK2": (6, 3), "IRF3": (6, 3),
       "ERK2": (6, 3), "p38": (6, 3), "Src": (5, -8), "AMPK": (6, -9),
       "PYGM": (-6, 7), "CDK1": (2, -11), "AKT1": (6, 3), "GSK3B": (6, 3),
       "PAK1": (6, 3), "TBK1": (6, 3), "IKBb": (6, 3)}
for r in rows:
    dx, dy = OFF.get(r["short"], (5, 3))
    ax1.annotate(r["short"], (r["mean_disp_exp"], r["mean_disp_pred"]),
                 textcoords="offset points", xytext=(dx, dy), fontsize=6.0)

ax1.set_xlim(0, lim); ax1.set_ylim(0, lim)
ax1.set_xlabel("Experimental displacement (Å)")
ax1.set_ylabel("Predicted displacement (Å)")
ax1.set_title("a   Magnitude of PTM-induced shift", loc="left", fontweight="bold")
ax1.legend(frameon=False, loc="upper left", fontsize=5.8, bbox_to_anchor=(0.0, 1.0))
ax1.spines[["top", "right"]].set_visible(False)
ax1.text(0.55, 0.97, "below 1:1 → under-estimate", transform=ax1.transAxes,
         fontsize=6, color="0.35", va="top", ha="center")

# ---------------- Panel b：方向余弦 / 幅度 r（AF3，可选 ESMFold） ----------------
rows_s = sorted(rows, key=lambda r: (r.get("mean_cosine") is None, r.get("mean_cosine") or 0))
labels = [r["short"] for r in rows_s]
yi = np.arange(len(rows_s))
ax2.barh(yi + 0.19, [r.get("mean_cosine") or 0 for r in rows_s], 0.36, color="#2ca02c",
         edgecolor="black", linewidth=0.4, label="Direction cosine (AF3)")
ax2.barh(yi - 0.19, [r.get("pearson_magnitude") or 0 for r in rows_s], 0.36, color="#ff7f0e",
         edgecolor="black", linewidth=0.4, label="Magnitude $r$ (AF3)")

# ESMFold 叠加（若可用）：在同一蛋白行上以空心菱形标注 cos / r，
# 与 AF3 条形对齐，避免 4 组条形过密。
# ESMFold overlay (if available): hollow diamonds for cos / r on the same protein
# rows, aligned with the AF3 bars (avoids an over-dense 4-bar grouping).
if esm:
    esm_by = {r["short"]: r for r in esm}
    ex_cos = [esm_by[l].get("mean_cosine") if l in esm_by else None for l in labels]
    ex_r = [esm_by[l].get("pearson_magnitude") if l in esm_by else None for l in labels]
    ax2.scatter([np.nan if v is None else v for v in ex_cos], yi + 0.19, marker="D", s=14,
                facecolor="none", edgecolor="#14599b", linewidth=0.8,
                label="Direction cosine (ESMFold)", zorder=5)
    ax2.scatter([np.nan if v is None else v for v in ex_r], yi - 0.19, marker="D", s=14,
                facecolor="none", edgecolor="#b35900", linewidth=0.8,
                label="Magnitude $r$ (ESMFold)", zorder=5)

ax2.axvline(0, color="0.5", lw=0.6)
ax2.set_yticks(yi); ax2.set_yticklabels(labels)
ax2.set_xlim(-0.45, 1.10)
ax2.set_xlabel("Agreement with experiment")
ax2.set_title("b   Direction vs magnitude agreement", loc="left", fontweight="bold")
ax2.legend(frameon=False, loc="lower right", ncol=1)
ax2.spines[["top", "right"]].set_visible(False)

fig.tight_layout(pad=0.8, w_pad=2.0)
fig.savefig(FIG / "Fig2_af3_vs_pdb.pdf", bbox_inches="tight")
fig.savefig(FIG / "Fig2_af3_vs_pdb.png", bbox_inches="tight", dpi=300)
plt.close(fig)
print(f"[OK] {FIG}/Fig2_af3_vs_pdb.pdf  (AF3 n={len(rows)}"
      + (f", ESMFold n={len(esm)}" if esm else ", ESMFold pending") + ")")
for r in rows:
    print(f"  {r['short']:8s} exp={r['mean_disp_exp']:.2f} pred={r['mean_disp_pred']:.2f} "
          f"cos={r.get('mean_cosine')} r={r.get('pearson_magnitude')}")
