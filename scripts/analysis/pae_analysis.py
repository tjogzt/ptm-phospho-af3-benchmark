#!/usr/bin/env python3
"""
pae_analysis.py — PAE（预测对齐误差）矩阵分析
===============================================
动机：目前只用了 pLDDT（单残基置信度）。但"磷酸化是否引起结构域相对移动"
本质上是**残基对之间**的问题，对应 AF3 的 PAE 矩阵，此前完全未使用。

指标：
  pae_global        全矩阵均值（整体结构确定性）
  pae_site_to_rest  位点区(±10 残基) → 其余残基 的平均 PAE（位点相对整体的取向置信度）
  pae_site_internal 位点区内部的平均 PAE
  delta_pae_*       phos 相对 unphos 的变化（正值 = 磷酸化后更不确定）

输出：results/analysis/pae_metrics.json + Fig5

用法：python3 scripts/analysis/pae_analysis.py
"""
import sys, json, zipfile
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

DOWNLOADS = _ROOT/"results"/"af3"/"downloads"
OUT = _ROOT/"results"/"analysis"
FIG = _ROOT/"results"/"figures"

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"axes.labelsize":7.5,
 "axes.titlesize":8,"xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":6.3,
 "axes.linewidth":0.6,"pdf.fonttype":42})

SITE_WIN = 10          # 位点窗口 ±10 残基


def collect_pae():
    """从全部 zip 中提取每作业的 PAE 指标（取 model_0，按作业聚合）"""
    recs = defaultdict(list)
    for zf in sorted(DOWNLOADS.glob("folds_*.zip")):
        if zf.name.startswith("._"): continue
        try:
            with zipfile.ZipFile(zf) as z:
                for name in z.namelist():
                    if "full_data" not in name: continue
                    if not name.endswith("_0.json"):   # 只取 model_0
                        continue
                    job = name.split("/")[0]
                    try:
                        d = json.loads(z.read(name))
                    except Exception:
                        continue
                    pae = d.get("pae")
                    resids = d.get("token_res_ids")
                    if not pae or not resids: continue
                    recs[job].append({
                        "pae": np.array(pae, dtype=float),
                        "resids": list(resids),
                        "plddt": d.get("atom_plddts"),
                    })
        except Exception:
            continue
    return recs


def pae_metrics(pae, resids, sites):
    """计算单模型的 PAE 指标"""
    n = len(resids)
    if pae.shape[0] != n: return None
    idx_site = [i for i, r in enumerate(resids)
                if any(abs(r - s) <= SITE_WIN for s in sites)]
    idx_rest = [i for i in range(n) if i not in set(idx_site)]
    out = {"pae_global": float(pae.mean())}
    if idx_site and idx_rest:
        out["pae_site_to_rest"] = float(pae[np.ix_(idx_site, idx_rest)].mean())
        out["pae_rest_to_site"] = float(pae[np.ix_(idx_rest, idx_site)].mean())
    if len(idx_site) > 1:
        out["pae_site_internal"] = float(pae[np.ix_(idx_site, idx_site)].mean())
    return out


def main():
    recs = collect_pae()
    print(f"收集到 {len(recs)} 个作业的 PAE")

    # 按蛋白/状态聚合
    by_prot = defaultdict(lambda: {"phos": [], "unphos": []})
    for job, items in recs.items():
        if "_seed" not in job: continue
        gene = job.split("_seed")[0]
        state = "phos" if gene.endswith("_phos") else "unphos"
        g = gene[:-5] if gene.endswith("_phos") else gene[:-7]
        for it in items:
            by_prot[g][state].append(it)

    # 映射 gene -> benchmark key
    gene2key = {v["gene"].lower(): k for k, v in BENCHMARK_PROTEINS.items()}

    rows = []
    for g, d in sorted(by_prot.items()):
        key = gene2key.get(g)
        if not key: continue
        sites = [int("".join(c for c in s if c.isdigit()))
                 for s in BENCHMARK_PROTEINS[key].get("phos_sites", [])]
        rec = {"gene": g, "protein": key, "sites": sites}
        for state in ("phos", "unphos"):
            ms = []
            for it in d[state]:
                m = pae_metrics(it["pae"], it["resids"], sites)
                if m: ms.append(m)
            if ms:
                for k in ms[0]:
                    rec[f"{state}_{k}"] = round(float(np.mean([m[k] for m in ms if k in m])), 2)
                rec[f"{state}_n"] = len(ms)
        # 变化量
        for k in ("pae_global", "pae_site_to_rest", "pae_site_internal"):
            a, b = rec.get(f"phos_{k}"), rec.get(f"unphos_{k}")
            if a is not None and b is not None:
                rec[f"delta_{k}"] = round(a - b, 2)
        rows.append(rec)

    # 合并性能指标
    comp = json.load(open(OUT/"af3_vs_pdb_summary.json"))
    perf = {}
    S2K = {"CDK2":"CDK2","ERK2":"ERK2","AURKA":"Aurora_A","IRF3":"IRF3","p38":"P38alpha",
           "CDK7":"CDK7","Src":"c-Src_kinase","AMPK":"AMPK_alpha2","PYGM":"Glycogen_Phosphorylase",
           "CDK1":"CDK1","ERK1":"ERK1","GSK3B":"GSK3B","PAK1":"PAK1","AKT1":"AKT1",
           "TBK1":"TBK1","IKBb":"IKBKB"}
    for r in comp["per_protein"]:
        k = S2K.get(r["short"])
        if k: perf[k] = r
    for rec in rows:
        p = perf.get(rec["protein"])
        if p:
            rec["amplitude_ratio"] = round(p["mean_disp_pred"]/p["mean_disp_exp"], 3) if p["mean_disp_exp"] else None
            rec["direction_cosine"] = p.get("mean_cosine")
            rec["magnitude_r"] = p.get("pearson_magnitude")

    (OUT/"pae_metrics.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False))

    print(f"\n{'蛋白':<8}{'pae全局(u)':>12}{'位点→其余(u)':>14}{'Δpae全局':>10}"
          f"{'Δ位点→其余':>12}{'幅度比':>9}")
    print("-"*76)
    for r in sorted(rows, key=lambda x: (x.get("amplitude_ratio") is None, x.get("amplitude_ratio") or 0)):
        def fv(k, d=0):
            v = r.get(k)
            return f"{v:>{d}.2f}" if v is not None else f"{'n/a':>{d}}"
        print(f"{r['gene']:<8}{fv('unphos_pae_global',12)}{fv('unphos_pae_site_to_rest',14)}"
              f"{fv('delta_pae_global',10)}{fv('delta_pae_site_to_rest',12)}"
              f"{fv('amplitude_ratio',9)}")

    # 相关分析
    valid = [r for r in rows if r.get("amplitude_ratio") is not None
             and r.get("unphos_pae_site_to_rest") is not None]
    if len(valid) >= 4:
        print(f"\n相关分析（n={len(valid)}）：")
        for feat, lab in [("unphos_pae_global","全局 PAE"),
                          ("unphos_pae_site_to_rest","位点→其余 PAE"),
                          ("unphos_pae_site_internal","位点内部 PAE"),
                          ("delta_pae_global","Δ全局 PAE"),
                          ("delta_pae_site_to_rest","Δ位点→其余 PAE")]:
            a = np.array([r[feat] for r in valid if r.get(feat) is not None])
            b = np.array([r["amplitude_ratio"] for r in valid if r.get(feat) is not None])
            if len(a) >= 4 and a.std() > 1e-9 and b.std() > 1e-9:
                print(f"  幅度比 vs {lab:<18}: r = {np.corrcoef(a,b)[0,1]:+.3f}")
    print(f"\n[OK] {OUT}/pae_metrics.json")
    return rows


def plot_leverage(rows, out_pdf, out_png):
    """Fig5：ΔPAE 与幅度比 —— 展示杠杆点问题"""
    import matplotlib.pyplot as plt
    valid = [r for r in rows if r.get("amplitude_ratio") is not None
             and r.get("delta_pae_site_to_rest") is not None]
    if len(valid) < 4: return
    x = np.array([r["delta_pae_site_to_rest"] for r in valid])
    y = np.array([r["amplitude_ratio"] for r in valid])
    names = [r["gene"].upper() for r in valid]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9), sharey=True)

    def draw(ax, mask, title):
        xs, ys, ns = x[mask], y[mask], [n for n, m in zip(names, mask) if m]
        ax.scatter(xs, ys, s=44, c="#1f77b4", edgecolor="black", linewidth=0.5, zorder=3)
        for a, b, nm in zip(xs, ys, ns):
            off = (-30, 2) if nm in ("IRF3",) else (6, 3)
            ax.annotate(nm, (a, b), textcoords="offset points", xytext=off, fontsize=6)
        if len(xs) > 3 and xs.std() > 1e-9:
            k, c0 = np.polyfit(xs, ys, 1)
            xr = np.linspace(xs.min()-0.5, xs.max()+0.5, 50)
            ax.plot(xr, k*xr+c0, ls="--", lw=0.8, color="0.55")
            rr = np.corrcoef(xs, ys)[0, 1]
            ax.text(0.03, 0.95, f"r = {rr:+.3f}  (n={len(xs)})", transform=ax.transAxes,
                    fontsize=6.5, va="top")
        ax.axhline(1.0, color="0.8", lw=0.5, ls=":")
        ax.set_xlabel(r"$\Delta$PAE (site$\rightarrow$rest): phos $-$ unphos")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)

    draw(ax1, np.ones(len(x), bool), "a   All pairs (incl. confounded)")
    draw(ax2, np.array([n != "IRF3" for n in names]), "b   After removing IRF3")
    ax1.set_ylabel("Amplitude ratio (AF3 / experiment)")

    fig.tight_layout(pad=0.8, w_pad=1.2)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[OK] {out_pdf}")


if __name__ == "__main__":
    _rows = main()
    plot_leverage(_rows, FIG / "Fig5_pae_leverage.pdf", FIG / "Fig5_pae_leverage.png")

