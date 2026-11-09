#!/usr/bin/env python3
"""
residue_spectrum.py — 逐残基位移谱分析：解释 AF3 为何低估幅度
===============================================================
两个核心问题：
  Q1 预测力是否随"距磷酸化位点距离"衰减？（局部作用 vs 远端作用）
  Q2 AF3 是否系统性压缩大幅度运动？（幅度压缩假说）

输入：results/analysis/af3_vs_pdb_per_residue.csv（由 compare_af3_vs_pdb.py 生成）
输出：results/figures/Fig3_residue_spectrum.pdf/.png
      results/analysis/residue_spectrum.json

用法：python3 scripts/analysis/residue_spectrum.py
"""
import json, csv
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys

_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

FIG = _ROOT/"results"/"figures"; FIG.mkdir(parents=True, exist_ok=True)
OUT = _ROOT/"results"/"analysis"; OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"axes.labelsize":7.5,
 "axes.titlesize":8,"xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":6.3,
 "axes.linewidth":0.6,"pdf.fonttype":42})

SHORT2KEY = {"CDK2":"CDK2","ERK2":"ERK2","AURKA":"Aurora_A","IRF3":"IRF3","p38":"P38alpha",
             "CDK7":"CDK7","Src":"c-Src_kinase","AMPK":"AMPK_alpha2",
             "PYGM":"Glycogen_Phosphorylase","CDK1":"CDK1"}


def main():
    rows = list(csv.DictReader(open(OUT/"af3_vs_pdb_per_residue.csv")))
    # 位点表
    sites = {}
    for k, v in BENCHMARK_PROTEINS.items():
        if not v.get("paired"): continue
        sites[v["gene"]] = [int("".join(c for c in s if c.isdigit()))
                            for s in v["phos_sites"]]

    # 补 distance-to-site
    recs = []
    for r in rows:
        key = SHORT2KEY.get(r["short"])
        if not key: continue
        gene = BENCHMARK_PROTEINS[key]["gene"]
        sl = sites.get(gene, [])
        if not sl: continue
        resid = int(r["resid"])
        d = min(abs(resid - s) for s in sl)
        recs.append({"short": r["short"], "resid": resid, "d_site": d,
                     "exp": float(r["disp_exp"]), "pred": float(r["disp_pred"])})

    print(f"有效残基 {len(recs)} 个（{len({x['short'] for x in recs})} 蛋白）")

    # ---------- Q1: 按距位点距离分箱 ----------
    bins = [(0,5),(5,10),(10,20),(20,40),(40,80),(80,10**6)]
    labels = ["0–5","5–10","10–20","20–40","40–80",">80"]
    q1 = []
    for (lo,hi), lab in zip(bins, labels):
        sub = [x for x in recs if lo <= x["d_site"] < hi]
        if not sub: continue
        e = np.array([x["exp"] for x in sub]); p = np.array([x["pred"] for x in sub])
        ratio = np.median(p)/np.median(e) if np.median(e) > 1e-6 else np.nan
        # 幅度相关
        corr = np.corrcoef(e, p)[0,1] if e.std()>1e-9 and p.std()>1e-9 else np.nan
        q1.append({"bin": lab, "n": len(sub), "median_exp": round(float(np.median(e)),3),
                   "median_pred": round(float(np.median(p)),3),
                   "ratio": round(float(ratio),3) if not np.isnan(ratio) else None,
                   "corr": round(float(corr),3) if not np.isnan(corr) else None})

    print("\nQ1  距磷酸化位点距离（残基）  n     exp中位  pred中位  比值  相关")
    print("-"*74)
    for x in q1:
        print(f"    {x['bin']:<12}{x['n']:>8}  {x['median_exp']:>7.2f}  "
              f"{x['median_pred']:>8.2f}  {x['ratio'] if x['ratio'] is not None else float('nan'):>6.2f}  "
              f"{x['corr'] if x['corr'] is not None else float('nan'):>6.3f}")

    # ---------- Q2: 按实验位移大小分箱（检验幅度压缩） ----------
    edges = [0,0.5,1,2,3,5,10**6]
    labs = ["<0.5","0.5–1","1–2","2–3","3–5",">5"]
    q2 = []
    for (lo,hi), lab in zip(zip(edges[:-1], edges[1:]), labs):
        sub = [x for x in recs if lo <= x["exp"] < hi]
        if not sub: continue
        e = np.array([x["exp"] for x in sub]); p = np.array([x["pred"] for x in sub])
        q2.append({"bin": lab, "n": len(sub), "median_exp": round(float(np.median(e)),3),
                   "median_pred": round(float(np.median(p)),3),
                   "ratio": round(float(np.median(p)/np.median(e)),3)})

    print("\nQ2  实验位移(Å)       n     exp中位  pred中位  比值")
    print("-"*56)
    for x in q2:
        print(f"    {x['bin']:<12}{x['n']:>8}  {x['median_exp']:>7.2f}  "
              f"{x['median_pred']:>8.2f}  {x['ratio']:>6.2f}")

    # ---------- 绘图 ----------
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 2.9))

    x1 = np.arange(len(q1))
    axA.bar(x1-0.19, [x["median_exp"] for x in q1], 0.36, color="#1f77b4",
            edgecolor="black", linewidth=0.4, label="Experiment")
    axA.bar(x1+0.19, [x["median_pred"] for x in q1], 0.36, color="#ff7f0e",
            edgecolor="black", linewidth=0.4, label="AF3")
    axA.set_xticks(x1); axA.set_xticklabels([x["bin"] for x in q1], fontsize=6.5)
    axA.set_xlabel("Distance to phosphosite (residues)")
    axA.set_ylabel("Median displacement (Å)")
    axA.set_title("a   Predictive power vs. distance from site", loc="left", fontweight="bold")
    axA.legend(frameon=False, loc="upper right")
    axA.spines[["top","right"]].set_visible(False)

    x2 = np.arange(len(q2))
    axB.plot(x2, [x["median_exp"] for x in q2], "o-", color="#1f77b4", ms=4,
             lw=1.0, label="Experiment")
    axB.plot(x2, [x["median_pred"] for x in q2], "s-", color="#ff7f0e", ms=4,
             lw=1.0, label="AF3")
    axB.plot([0, len(q2)-1], [0, q2[-1]["median_exp"]], ls=":", lw=0.7, color="0.6")
    axB.set_xticks(x2); axB.set_xticklabels([x["bin"] for x in q2], fontsize=6.5)
    axB.set_xlabel("Experimental displacement bin (Å)")
    axB.set_ylabel("Median displacement (Å)")
    axB.set_title("b   Amplitude compression", loc="left", fontweight="bold")
    axB.legend(frameon=False, loc="upper left")
    axB.spines[["top","right"]].set_visible(False)
    # 标注压缩比
    for i, x in enumerate(q2):
        if x["ratio"] < 0.5:
            axB.annotate(f'{x["ratio"]:.2f}×', (i, x["median_pred"]),
                         textcoords="offset points", xytext=(0, -9),
                         fontsize=5.5, ha="center", color="#ff7f0e")

    fig.tight_layout(pad=0.8, w_pad=1.8)
    fig.savefig(FIG/"Fig3_residue_spectrum.pdf", bbox_inches="tight")
    fig.savefig(FIG/"Fig3_residue_spectrum.png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    summary = {"n_residues": len(recs), "by_distance": q1, "by_amplitude": q2,
               "overall": {
                   "median_exp": round(float(np.median([x["exp"] for x in recs])),3),
                   "median_pred": round(float(np.median([x["pred"] for x in recs])),3),
                   "ratio": round(float(np.median([x["pred"] for x in recs]) /
                                        max(np.median([x["exp"] for x in recs]),1e-9)),3),
               }}
    (OUT/"residue_spectrum.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"\n[OK] {FIG}/Fig3_residue_spectrum.pdf")
    print(f"[OK] {OUT}/residue_spectrum.json")
    return summary


if __name__ == "__main__":
    main()
