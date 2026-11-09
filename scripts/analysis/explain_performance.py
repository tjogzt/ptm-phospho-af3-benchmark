#!/usr/bin/env python3
"""
explain_performance.py — 解释 AF3 预测表现差异的成因
======================================================
Why does AF3 perform better on some phospho-pairs than others?

检验两类假设 / Two hypotheses:
H1 位移聚集性：实验位移若集中在少数残基（局部重排），AF3 更易捕捉；
               若分散于整个结构域（整体刚体运动），单体预测难以复现。
   指标：top-10% 残基贡献的位移能量占比（集中度高 → 局部）
H2 位点区柔性：位点周围若本就低置信度（柔性 loop），AF3 的响应会更受限。
   指标：AF3 unphos 在位点 ±SITE_WIN 残基的平均 pLDDT

【v3 (B2/B6) 修复】
  - B2：pLDDT 提取改用氨基酸类残基判定（含以 HETATM 出现的 SEP/TPO/PTR），
        与 compare/geometry 的 Cα 提取口径一致；位点窗统一 ±10。
  - B6：分析队列改为冻结 manifest 的 primary cohort（不再用旧的
        pair_audit.json PASS 判定），逐蛋白指标按 benchmark key 从
        af3_vs_pdb_summary.json 读取，避免标签错配。
  - 清理死代码：原 plot() 定义在 `if __name__ == "__main__"` 之后从未被调用；
    现纳入 main() 调用链，生成 Fig4。

历史（供追溯）/ History:
  v1 (2026-10-25): 6 组纯净配对；位点 pLDDT vs 幅度比 r=+0.649（未显著）
  v2 (2026-10-27): 扩样至 n=11；r 降至 +0.310 → 原假设未获支持，已撤回
    详见 audit/FIG4_HYPOTHESIS_REFUTED.md
    表述要求：写作时应称为"观察到的趋势/机制假设"，不可表述为确立的相关性
  v3 (2026-11-08): B2/B6 重算；队列与标签与 manifest/compare 统一

输出：results/analysis/performance_explained.json + results/figures/Fig4_*.pdf
用法：python3 scripts/analysis/explain_performance.py
"""
import sys, json, csv
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser

_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_dataset import BENCHMARK_PROTEINS
from struct_utils import aa_from_resname        # B2：氨基酸类残基判定（含 HETATM）

CIF = _ROOT / "results" / "af3" / "models" / "v2"
OUT = _ROOT / "results" / "analysis"
FIG = _ROOT / "results" / "figures"

SITE_WIN = 10        # 位点窗 ±10（B2：与 compare/geometry 统一）

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 7.5,
                     "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 6.3, "axes.linewidth": 0.6, "pdf.fonttype": 42})

MANIFEST = _ROOT / "results" / "manifest_frozen.json"
if MANIFEST.exists():
    _mf = json.loads(MANIFEST.read_text())
    PAIRED = list(_mf["cohorts"]["primary"])
    PURE = set(_mf["cohorts"]["pure"])
else:
    PAIRED = [k for k, v in BENCHMARK_PROTEINS.items() if v.get("paired")]
    PURE = set()

SITE = {k: [int("".join(c for c in s if c.isdigit())) for s in v.get("phos_sites", [])]
        for k, v in BENCHMARK_PROTEINS.items()}


def plddt_profile(cif_path):
    """从 AF3 CIF 提取 per-residue pLDDT（B-factor 列）。

    纳入以 HETATM 出现的修饰残基（SEP/TPO/PTR），与 Cα 提取口径一致（B2）。
    Includes modified residues recorded as HETATM (SEP/TPO/PTR), matching the
    Cα-extraction scope used elsewhere (B2).
    """
    p = MMCIFParser(QUIET=True)
    try:
        s = p.get_structure("x", str(cif_path))
    except Exception:
        return {}
    model = next(iter(s))
    out = {}
    for ch in model:
        for res in ch:
            if res.id[0][0] not in (" ", "H"):
                continue
            if aa_from_resname(res.get_resname()) is None:
                continue
            vals = [a.get_bfactor() for a in res]
            if vals:
                out[res.id[1]] = float(np.mean(vals))
        if out:
            break
    return out


def main():
    comp = json.loads((OUT / "af3_vs_pdb_summary.json").read_text())
    perf = {r["protein"]: r for r in comp["per_protein"]
            if r.get("mean_disp_exp") and r.get("mean_disp_pred") is not None}

    rows = list(csv.DictReader(open(OUT / "af3_vs_pdb_per_residue.csv")))
    by_prot = defaultdict(list)
    for r in rows:
        by_prot[r["protein"]].append(r)

    results = []
    for key in PAIRED:
        info = BENCHMARK_PROTEINS.get(key)
        if not info or not info.get("paired"):
            continue
        sites = SITE.get(key, [])
        recs = by_prot.get(key, [])
        if not recs:
            continue
        exp = np.array([float(r["disp_exp"]) for r in recs])
        if exp.sum() <= 0:
            continue

        # H1 聚集性：按位移降序，top-10% 残基贡献的能量占比
        e2 = exp ** 2
        k = max(1, int(np.ceil(0.10 * len(exp))))
        top_frac = float(np.sort(e2)[::-1][:k].sum() / e2.sum())

        # H2 位点区 pLDDT（AF3 unphos，用 model 0 seed1）
        gene = info["gene"].lower()
        cif_u = sorted((CIF / gene / "unphos").glob("*seed1_model_0.cif"))
        cif_p = sorted((CIF / gene / "phos").glob("*seed1_model_0.cif"))
        prof_u = plddt_profile(cif_u[0]) if cif_u else {}
        prof_p = plddt_profile(cif_p[0]) if cif_p else {}

        def win(prof):
            if not prof or not sites:
                return None
            v = [prof.get(r) for r in range(min(sites) - SITE_WIN, max(sites) + SITE_WIN + 1)
                 if prof.get(r) is not None]
            return float(np.mean(v)) if v else None

        site_plddt_u, site_plddt_p = win(prof_u), win(prof_p)
        global_u = float(np.mean(list(prof_u.values()))) if prof_u else None

        r = perf.get(key)
        if r is None:
            continue
        results.append({
            "protein": key, "short": r["short"],
            "effect_type": next((p.get("effect_type") for p in _mf.get("proteins", [])
                                 if p["key"] == key), "n/a") if MANIFEST.exists() else "n/a",
            "site_plddt_unphos": round(site_plddt_u, 1) if site_plddt_u else None,
            "site_plddt_phos": round(site_plddt_p, 1) if site_plddt_p else None,
            "global_plddt_unphos": round(global_u, 1) if global_u else None,
            "site_minus_global": round(site_plddt_u - global_u, 1)
                                 if (site_plddt_u and global_u) else None,
            "disp_top10pct_energy_frac": round(top_frac, 3),
            "amplitude_ratio": round(r["mean_disp_pred"] / r["mean_disp_exp"], 3),
            "direction_cosine": r.get("mean_cosine"),
            "magnitude_r": r.get("pearson_magnitude"),
        })

    results.sort(key=lambda x: x["amplitude_ratio"])
    print("AF3 表现差异的成因分析")
    print("=" * 104)
    print(f"{'蛋白':<8}{'位点pLDDT(u)':>13}{'全局pLDDT':>11}{'差值':>8}"
          f"{'位移集中度':>11}{'幅度比':>9}{'方向cos':>9}")
    print("-" * 104)
    for x in results:
        print(f"{x['short']:<8}"
              f"{(x['site_plddt_unphos'] if x['site_plddt_unphos'] else 0):>13.1f}"
              f"{(x['global_plddt_unphos'] if x['global_plddt_unphos'] else 0):>11.1f}"
              f"{(x['site_minus_global'] if x['site_minus_global'] is not None else 0):>8.1f}"
              f"{x['disp_top10pct_energy_frac']:>11.2f}"
              f"{(x['amplitude_ratio'] or 0):>9.2f}"
              f"{(x['direction_cosine'] or 0):>9.3f}")

    # 相关分析（探索性；须在写作中表述为"观察到的趋势"）
    vals = [x for x in results if x["amplitude_ratio"] is not None
            and x["site_plddt_unphos"] is not None]
    print(f"\n相关分析（n={len(vals)}，探索性；不得表述为确立的相关性）：")
    if len(vals) >= 4:
        for feat, lab in [("site_plddt_unphos", "位点pLDDT"),
                          ("disp_top10pct_energy_frac", "位移集中度"),
                          ("site_minus_global", "位点-全局pLDDT差")]:
            sub = [x for x in vals if x.get(feat) is not None]
            a = np.array([x[feat] for x in sub])
            b = np.array([x["amplitude_ratio"] for x in sub])
            if a.std() > 1e-9 and b.std() > 1e-9:
                rr = float(np.corrcoef(a, b)[0, 1])
                print(f"  幅度比 vs {lab}: r = {rr:+.3f} (n={len(sub)})")

    (OUT / "performance_explained.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))
    print(f"\n[OK] {OUT}/performance_explained.json")

    try:
        plot(results, FIG / "Fig4_plddt_vs_amplitude.pdf",
             FIG / "Fig4_plddt_vs_amplitude.png")
    except Exception as e:
        print(f"[WARN] Fig4 生成失败: {type(e).__name__}: {e}")
    return results


def plot(results, out_pdf, out_png):
    """Fig4：位点区 AF3 置信度 vs 幅度预测能力（探索性，标题中性）。

    Fig4 — site-region AF3 confidence vs amplitude accuracy. Deliberately
    descriptive titles; the relationship is exploratory and was previously
    refuted on a smaller cohort.
    """
    vals = [x for x in results
            if x["amplitude_ratio"] is not None and x["site_plddt_unphos"] is not None]
    if len(vals) < 3:
        return
    x = np.array([v["site_plddt_unphos"] for v in vals])
    y = np.array([v["amplitude_ratio"] for v in vals])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9))

    ax1.scatter(x, y, s=46, c="#1f77b4", edgecolor="black", linewidth=0.5, zorder=3)
    OFF = {"CDK2": (-24, -2), "p38": (6, 2), "ERK2": (6, 2),
           "PYGM": (6, 2), "Src": (-22, 1), "AURKA": (7, -1)}
    for v in vals:
        dx, dy = OFF.get(v["short"], (5, 3))
        ax1.annotate(v["short"], (v["site_plddt_unphos"], v["amplitude_ratio"]),
                     textcoords="offset points", xytext=(dx, dy), fontsize=6.3)
    if x.std() > 1e-9:
        k, b = np.polyfit(x, y, 1)
        xr = np.linspace(x.min() - 1.5, x.max() + 1.5, 50)
        ax1.plot(xr, k * xr + b, ls="--", lw=0.8, color="0.55")
        r = float(np.corrcoef(x, y)[0, 1])
        ax1.text(0.03, 0.95, f"r = {r:+.2f}  (n={len(vals)}, exploratory)",
                 transform=ax1.transAxes, fontsize=6.5, va="top")
    ax1.axhline(1.0, color="0.75", lw=0.5, ls=":")
    ax1.set_xlabel("AF3 pLDDT at phosphosite region")
    ax1.set_ylabel("Amplitude ratio (AF3 / experiment)")
    ax1.set_title("a   Site-region confidence vs amplitude ratio (exploratory)",
                  loc="left", fontweight="bold")
    ax1.spines[["top", "right"]].set_visible(False)

    # Panel b：位点 vs 全局 pLDDT
    labs = [v["short"] for v in vals]
    sy = np.arange(len(vals))
    ax2.barh(sy + 0.19, [v["global_plddt_unphos"] for v in vals], 0.36,
             color="#7f7f7f", edgecolor="black", linewidth=0.4, label="Global")
    ax2.barh(sy - 0.19, [v["site_plddt_unphos"] for v in vals], 0.36,
             color="#d62728", edgecolor="black", linewidth=0.4, label="Site region")
    ax2.set_yticks(sy)
    ax2.set_yticklabels(labs)
    ax2.set_xlabel("AF3 pLDDT (unphos)")
    ax2.set_title("b   Phosphosite regions tend to be low-confidence",
                  loc="left", fontweight="bold")
    ax2.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
    ax2.set_xlim(0, 105)
    ax2.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(pad=0.8, w_pad=1.8)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[OK] {out_pdf}")


if __name__ == "__main__":
    main()
