#!/usr/bin/env python3
"""
geometry_decomposition.py — 位移几何分解（v3，分块叠合方案 + B2/B3 修复）
======================================================================
Geometry decomposition of phosphorylation-induced displacement (v3).

【方案 / Method】
把两态共有残基按序列切成 k 块，每块各自叠合到参考结构：
  r_full    全残基一次叠合后的 RMSD        —— 整体形变
  r_block   k 块各自叠合后的块内 RMSD 均值 —— 局部柔性
  hinge     = 1 - r_block / r_full          —— 块间相对运动（铰链）占比
              hinge 高 → 差异主要来自块与块的相对位置变化（铰链/域运动）
              hinge 低 → 差异主要是块内形变（局部重排/柔性）
  disp_cv   位移大小的变异系数（离散度）

  检验假设：AF3 若擅长局部形变、不擅长铰链运动，则 hinge 应与幅度比负相关。
  注：当前证据不足（不显著），仅作探索性描述，不得表述为确立的相关性。

【v1 的错误】在全残基 Kabsch 叠合后计算方向一致性 |Σd|/Σ|d|：叠合本身令质心
对齐，故 Σd ≡ 0，该指标恒为零（数学上必然）。

【v3 (B2/B3) 修复】
  - B2：Cα 提取改用 struct_utils.load_ca，纳入以 HETATM 出现的修饰残基
        （SEP/TPO/PTR）。旧代码用 `res.id[0] != " "` 会静默丢弃磷酸化位点
        残基，使度量在两侧都缺失位点。位点窗统一为 ±SITE_WIN(=10)。
  - B3：链选择与 compare_af3_vs_pdb.py 同源（pdb_chains 逐结构解析、缺失时
        回退到首链），并记录回退事件到 chain_resolution_log，消除同一蛋白在
        两条管线上的数值矛盾。
  - 清理死代码：原 plot_core 定义在 `if __name__ == "__main__"` 之后从未被
    调用；现已纳入 main() 调用链。Fig6 标题改为中性、可辩护表述。

输出：results/analysis/geometry_metrics.json + results/figures/Fig6_*.pdf
用法：python3 scripts/analysis/geometry_decomposition.py
"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from difflib import SequenceMatcher

_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_dataset import BENCHMARK_PROTEINS
from struct_utils import load_ca                 # B2/B3：共享 Cα 提取 + 链选择

PDBDIR = _ROOT / "data" / "pdb"
OUT = _ROOT / "results" / "analysis"
FIG = _ROOT / "results" / "figures"

SITE_WIN = 10        # 位点窗 ±10 残基（B2：与 compare / pae / explain 统一）
N_BLOCKS = 3         # 分块叠合的块数

# 分析队列来自冻结 manifest（B1），保证与论文 Table 1 同源
MANIFEST = _ROOT / "results" / "manifest_frozen.json"
if MANIFEST.exists():
    _mf = json.loads(MANIFEST.read_text())
    PAIRED = list(_mf["cohorts"]["primary"])
    PURE = set(_mf["cohorts"]["pure"])
    EFFECT = {p["key"]: p.get("effect_type", "n/a") for p in _mf.get("proteins", [])}
else:                                                       # 回退
    PAIRED = [k for k, v in BENCHMARK_PROTEINS.items() if v.get("paired")]
    PURE, EFFECT = set(), {}

SITE = {k: [int("".join(c for c in s if c.isdigit())) for s in v.get("phos_sites", [])]
        for k, v in BENCHMARK_PROTEINS.items()}

CHAIN_LOG = []       # B3：链回退审计

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 7.5,
                     "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 6.3, "axes.linewidth": 0.6, "pdf.fonttype": 42})


def ca_chain(path, requested):
    """读 Cα（纳入 SEP/TPO/PTR）；链选择与 compare 同源并记录回退（B2/B3）。

    Read Cα (including SEP/TPO/PTR); chain selection mirrors the compare
    pipeline and records any fallback (B2/B3).
    """
    coords, seq, resids, meta = load_ca(path, chain=requested)
    CHAIN_LOG.append({"file": Path(path).name, "requested": requested, **meta})
    return coords, seq, resids


def kabsch(P, Q):
    """求 (R, t) 使 P @ R + t ≈ Q（行向量约定）。

    本机 numpy 以 Apple Accelerate 作 BLAS，其 (3×N)·(N×3) 核会对有限输入
    伪设 divide/overflow/underflow FPU 标志并报告为 RuntimeWarning（数值与
    纯 numpy einsum 完全一致，已核实）。故在 errstate 中屏蔽，不影响结果。
    """
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        H = Pc.T @ Qc
        U, _, Vt = np.linalg.svd(H)
        det = float(np.linalg.det(U @ Vt))
        D = np.diag([1.0, 1.0, 1.0 if det >= 0 else -1.0])
        R = U @ D @ Vt
        t = Q.mean(0) - P.mean(0) @ R
    return R, t


def rmsd_after_fit(A, B):
    """A 叠合到 B 后的 RMSD（<3 残基返回 None）。"""
    if len(A) < 3:
        return None
    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        R, t = kabsch(A, B)
        return float(np.sqrt((((A @ R + t) - B) ** 2).sum(1).mean()))


def align_pairs(a, b):
    """序列对齐 → [(idxA, idxB)] 匹配对。"""
    sm = SequenceMatcher(None, a, b, autojunk=False)
    return [(blk.a + k, blk.b + k)
            for blk in sm.get_matching_blocks() for k in range(blk.size)]


def resolve_chains(info, up_id, ph_id):
    """逐结构解析链（与 compare 同源；CDK7: 7B5O 在链 J、8PYR 在链 A）。"""
    chs = info.get("pdb_chains", {})
    chain = info.get("chain", "A")
    ch_up = chs.get(up_id, chs.get(up_id.upper(), chain))
    ch_ph = chs.get(ph_id, chs.get(ph_id.upper(), chain))
    return ch_up, ch_ph


def decompose(key, info):
    """对单个蛋白做几何分解；返回 dict 或 None（数据不足）。"""
    pdbs = info["pdbs"]
    up_id = next((x for x in pdbs.get("unphos", [])
                  if (PDBDIR / f"{x.upper()}.pdb").exists()), None)
    ph_id = next((x for x in pdbs.get("phos", [])
                  if (PDBDIR / f"{x.upper()}.pdb").exists()), None)
    if not (up_id and ph_id):
        return None
    ch_up, ch_ph = resolve_chains(info, up_id, ph_id)
    C_up, s_up, id_up = ca_chain(PDBDIR / f"{up_id.upper()}.pdb", ch_up)
    C_ph, s_ph, id_ph = ca_chain(PDBDIR / f"{ph_id.upper()}.pdb", ch_ph)
    if len(C_up) < 60 or len(C_ph) < 60:
        return None
    m = align_pairs(s_up, s_ph)
    if len(m) < 60:
        return None
    ia = np.array([x[0] for x in m])
    ib = np.array([x[1] for x in m])
    P = C_up[ia]
    Q = C_ph[ib]
    resids = [id_up[i] for i in ia]

    # 整体叠合
    r_full = rmsd_after_fit(Q, P)
    if not r_full or r_full <= 0:
        return None

    # 分块叠合：按序列切 N_BLOCKS 段
    n = len(P)
    bounds = np.linspace(0, n, N_BLOCKS + 1).astype(int)
    block_r = []
    for i in range(N_BLOCKS):
        lo, hi = bounds[i], bounds[i + 1]
        if hi - lo < 10:
            continue
        r = rmsd_after_fit(Q[lo:hi], P[lo:hi])
        if r is not None:
            block_r.append(r)
    if not block_r:
        return None
    r_block = float(np.mean(block_r))
    hinge = 1.0 - (r_block / r_full) if r_full > 1e-9 else np.nan

    # 叠合后位移大小
    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        R, t = kabsch(Q, P)
        d = (Q @ R + t) - P
    mag = np.linalg.norm(d, axis=1)

    # 位点窗 vs 其余（±SITE_WIN）
    sites = SITE.get(key, [])
    site_disp = rest_disp = None
    if sites:
        near = [i for i, r in enumerate(resids) if any(abs(r - s) <= SITE_WIN for s in sites)]
        if len(near) >= 3:
            rest = [i for i in range(len(resids)) if i not in set(near)]
            site_disp = round(float(mag[near].mean()), 3)
            rest_disp = round(float(mag[rest].mean()), 3) if rest else None

    return {
        "protein": key, "gene": info["gene"], "short": _short(key, info),
        "effect_type": EFFECT.get(key, "n/a"),
        "pure_pair": (key in PURE),
        "unphos_pdb": up_id, "phos_pdb": ph_id,
        "unphos_chain": ch_up, "phos_chain": ch_ph,
        "n_res": int(n),
        "r_full": round(r_full, 3), "r_block": round(r_block, 3),
        "hinge": round(float(hinge), 3) if not np.isnan(hinge) else None,
        "mean_disp": round(float(mag.mean()), 3),
        "max_disp": round(float(mag.max()), 3),
        "disp_cv": round(float(mag.std() / mag.mean()), 3) if mag.mean() > 1e-9 else None,
        "site_mean_disp": site_disp, "rest_mean_disp": rest_disp,
    }


def _short(key, info):
    """compare 脚本中的简短标签（与 af3_vs_pdb_summary.json 的 SHORT 对齐）。"""
    return {"Aurora_A": "AURKA", "P38alpha": "p38", "c-Src_kinase": "Src",
            "AMPK_alpha2": "AMPK", "Glycogen_Phosphorylase": "PYGM",
            "IKBKB": "IKBb", "ERK2": "ERK2", "ERK1": "ERK1"}.get(key, info["gene"])


def _perf_map():
    """从 af3_vs_pdb_summary.json 读逐蛋白指标（键为 benchmark key）。"""
    f = OUT / "af3_vs_pdb_summary.json"
    if not f.exists():
        return {}
    comp = json.loads(f.read_text())
    return {r["protein"]: r for r in comp.get("per_protein", [])}


def main():
    perf = _perf_map()
    rows = []
    print(f"{'蛋白':<8}{'r_full':>8}{'r_block':>9}{'hinge':>8}{'位移CV':>8}"
          f"{'平均位移':>9}{'幅度比':>9}{'cos':>8}{'位点/其余':>11}")
    print("-" * 82)
    for key in PAIRED:
        info = BENCHMARK_PROTEINS.get(key)
        if not info or not info.get("paired"):
            continue
        try:
            g = decompose(key, info)
        except Exception as e:
            print(f"{key:<8} 失败: {type(e).__name__}: {e}")
            continue
        if not g:
            print(f"{key:<8} 数据不足，跳过")
            continue
        p = perf.get(key)
        if p:
            g["short"] = p.get("short", g["short"])      # 与 compare 摘要标签对齐
            g["mean_disp_exp"] = p.get("mean_disp_exp")
            g["mean_disp_pred"] = p.get("mean_disp_pred")
            e = p.get("mean_disp_exp")
            g["amplitude_ratio"] = (round(p["mean_disp_pred"] / e, 3)
                                    if e else None)
            g["direction_cosine"] = p.get("mean_cosine")
            g["pearson_magnitude"] = p.get("pearson_magnitude")
        rows.append(g)
        cos_s = f"{g.get('direction_cosine'):+.3f}" if g.get("direction_cosine") is not None else "  n/a"
        site_s = (f"{g['site_mean_disp']:.2f}/{g['rest_mean_disp']:.2f}"
                  if g.get("site_mean_disp") is not None and g.get("rest_mean_disp") is not None
                  else "   n/a")
        print(f"{g['short']:<8}{g['r_full']:>8.2f}{g['r_block']:>9.2f}"
              f"{(g['hinge'] if g['hinge'] is not None else 0):>8.3f}"
              f"{(g['disp_cv'] or 0):>8.2f}{g['mean_disp']:>9.2f}"
              f"{(g.get('amplitude_ratio') or 0):>9.2f}{cos_s:>8}{site_s:>11}")

    (OUT / "geometry_metrics.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False))

    valid = [r for r in rows if r.get("amplitude_ratio") is not None and r.get("hinge") is not None]
    print(f"\n相关分析（n={len(valid)}，探索性）：")
    for feat, lab in [("hinge", "hinge指数(铰链占比)"), ("r_full", "整体形变 RMSD"),
                      ("r_block", "块内柔性 RMSD"), ("disp_cv", "位移离散度"),
                      ("mean_disp", "平均位移")]:
        sub = [r for r in valid if r.get(feat) is not None]
        if len(sub) < 5:
            continue
        a = np.array([r[feat] for r in sub])
        b = np.array([r["amplitude_ratio"] for r in sub])
        if a.std() < 1e-9 or b.std() < 1e-9:
            continue
        rr = float(np.corrcoef(a, b)[0, 1])
        n = len(sub)
        t = rr * np.sqrt((n - 2) / max(1 - rr ** 2, 1e-9))
        try:
            from scipy import stats
            p = float(2 * stats.t.sf(abs(t), df=n - 2))
        except ImportError:
            from math import erf, sqrt
            p = float(2 * (1 - 0.5 * (1 + erf(abs(t) / sqrt(2)))))
        print(f"  幅度比 vs {lab:<18}: r = {rr:+.3f}  (n={n}, p={p:.3f})")
    print(f"\n[OK] {OUT}/geometry_metrics.json")

    try:
        plot_core(rows, FIG / "Fig6_geometry_decomposition.pdf",
                  FIG / "Fig6_geometry_decomposition.png")
    except Exception as e:
        print(f"[WARN] Fig6 生成失败: {type(e).__name__}: {e}")

    if CHAIN_LOG:      # B3：链回退审计
        (OUT / "chain_resolution_log_geometry.json").write_text(
            json.dumps(CHAIN_LOG, indent=2, ensure_ascii=False))
    return rows


def plot_core(rows, out_pdf, out_png):
    """Fig6：几何分解汇总 —— 幅度比/方向余弦（按效应类型分层）+ hinge 探索。

    Fig6 — geometry summary: amplitude ratio & direction cosine (stratified by
    effect type) plus the exploratory hinge relationship. Titles are deliberately
    descriptive to avoid over-claiming.
    """
    recs = [r for r in rows if r.get("amplitude_ratio") is not None]
    if len(recs) < 5:
        return
    recs.sort(key=lambda r: r["amplitude_ratio"])
    labs = [r["short"] for r in recs]
    ratio = np.array([r["amplitude_ratio"] for r in recs])
    cos = np.array([r.get("direction_cosine") if r.get("direction_cosine") is not None else np.nan
                    for r in recs])
    EFFCOL = {"local_rearrangement": "#1f77b4", "allosteric": "#ff7f0e",
              "disorder_order": "#2ca02c"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 3.2))

    # Panel a：幅度比（条形，按效应类型着色）+ 方向余弦（点）
    y = np.arange(len(recs))
    colors = [EFFCOL.get(r.get("effect_type"), "0.6") for r in recs]
    ax1.barh(y, ratio, height=0.62, color=colors, edgecolor="black", linewidth=0.4)
    ax1.axvline(1.0, color="#d62728", lw=0.9, ls="--")
    ax1.text(1.0, len(recs) - 0.2, " parity", color="#d62728", fontsize=6, va="top")
    ax1.scatter(cos, y, s=22, marker="D", facecolor="white",
                edgecolor="black", linewidth=0.6, zorder=5)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labs)
    ax1.set_xlabel("Amplitude ratio (bars) / direction cosine (diamonds)")
    ax1.set_title("a   AF3 under-shoots amplitude and shows mixed direction",
                  loc="left", fontweight="bold")
    ax1.spines[["top", "right"]].set_visible(False)

    handles = [plt.Line2D([0], [0], marker="s", ls="", color=c, label=k)
               for k, c in EFFCOL.items()]
    ax1.legend(handles=handles, frameon=False, loc="lower right", fontsize=5.6)

    # Panel b：hinge vs 幅度比（探索性）
    hinge = np.array([r["hinge"] for r in recs])
    ax2.scatter(hinge, ratio, s=40, c="#2ca02c", edgecolor="black", linewidth=0.5, zorder=3)
    for r in recs:
        ax2.annotate(r["short"], (r["hinge"], r["amplitude_ratio"]),
                     textcoords="offset points", xytext=(5, 3), fontsize=5.8)
    if hinge.std() > 1e-9:
        k2, c2 = np.polyfit(hinge, ratio, 1)
        xr = np.linspace(hinge.min() - 0.02, hinge.max() + 0.02, 50)
        ax2.plot(xr, k2 * xr + c2, ls="--", lw=0.8, color="0.55")
        r2 = float(np.corrcoef(hinge, ratio)[0, 1])
        ax2.text(0.05, 0.06, f"r = {r2:+.3f}  (n={len(recs)}, exploratory)",
                 transform=ax2.transAxes, fontsize=6.5, va="bottom")
    ax2.axhline(1.0, color="0.75", lw=0.5, ls=":")
    ax2.set_xlabel("Hinge index (inter-block / total motion)")
    ax2.set_ylabel("Amplitude ratio (AF3 / experiment)")
    ax2.set_title("b   Hinge-type motion vs amplitude ratio (exploratory)",
                  loc="left", fontweight="bold")
    ax2.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(pad=0.8, w_pad=1.6)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[OK] {out_pdf}")


if __name__ == "__main__":
    main()
