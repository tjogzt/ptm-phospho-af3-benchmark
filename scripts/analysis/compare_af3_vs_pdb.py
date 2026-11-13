#!/usr/bin/env python3
"""
compare_af3_vs_pdb.py — P1：AF3 预测位移 vs PDB 实验位移
==========================================================
核心问题：AF3 引入磷酸化后预测出的构象变化，是否与晶体学观察到的变化一致？

方法（统一坐标系下比较位移向量）：
  1. 以 PDB unphos 为参考框架 Q
  2. 将 PDB phos 叠合到 Q（共同 Cα，Kabsch）→ Δ_exp = PDB_phos_aligned − Q
  3. AF3 两态**各自独立**叠合到同一参考框架 Q（每次预测都有各自的任意坐标系，
     故 unphos、phos 分别 fit→Q）：
     Δ_pred = (AF3_phos fit→Q) − (AF3_unphos fit→Q)
     独立叠合消除整体刚体位移；若只把 unphos 的变换套用到 phos，两态预测的
     坐标系差异会被误当作磷酸化位移（曾导致 25–39 Å 假信号）。
  4. 在共同残基上比较 Δ_pred 与 Δ_exp：
       - 位移幅度的 Pearson 相关
       - 位移方向的余弦相似度（逐残基，取位移 > 1 Å 者）
       - top-N 位移残基的重合度

投影约定敏感性 / Projection sensitivity（W1）：
  `indep`  = 正文默认协议（两态各自独立叠合到参考框架，见步骤 3）
  `direct` = 两个 AF3 态先互相叠合（AF3 phos → AF3 unphos）再取差，
             度量模型自身的构象变化，不参照实验坐标。
  注：AF3 的两态是两个**独立预测**，各自处于任意坐标系；因此任何直接做
  (P_phos − P_unphos) 的"共享变换"约定在物理上无定义（差值混入了整体
  刚体位移，实测会给出 ~28 Å 的伪信号），故不纳入。可行的约定只有上述两种。
  两者下队列级预测位移与 r 一并写入 projection_sensitivity.json，
  用于证明「近常量预测位移」这一结论不依赖于叠合约定。

输出：
  results/figures/Fig2_af3_vs_pdb_displacement.pdf/.png
  results/analysis/af3_vs_pdb_per_residue.csv
  results/analysis/af3_vs_pdb_summary.json
  results/analysis/projection_sensitivity.json

用法：python scripts/analysis/compare_af3_vs_pdb.py [--model K]
  默认聚合每蛋白全部 seed×model（通常 15 个）；--model K 仅用第 K 号 model。
"""
import sys, json, argparse
from pathlib import Path
import numpy as np
from Bio.SeqUtils import seq1
from difflib import SequenceMatcher
from struct_utils import load_ca          # B2/B3：共享 Cα 提取（纳入 SEP/TPO/PTR）+ 链选择

_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

CIF = _ROOT / "results" / "af3" / "models" / "v2"   # 仅用经位点核准的 v2 数据
PDBDIR = _ROOT / "data" / "pdb"
OUT = _ROOT / "results" / "analysis"
FIG = _ROOT / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)

SITE_WIN = 10        # 位点窗 ±10 残基（B2：与 pae_analysis/explain_performance 统一）

# 分析队列来自**冻结 manifest**（B1），保证与论文 Table 1 同源。
# The analysis cohort comes from the frozen manifest (B1) so it is
# identical to the one reported in the paper (Table 1).
MANIFEST = _ROOT / "results" / "manifest_frozen.json"
if MANIFEST.exists():
    _mf = json.loads(MANIFEST.read_text())
    PAIRED = list(_mf["cohorts"]["primary"])          # paired + 两态 AF3 齐全
    PURE = set(_mf["cohorts"]["pure"])                # P1–P5 全过
else:                                                 # 回退：直接从 benchmark 推导
    PAIRED = [k for k, v in BENCHMARK_PROTEINS.items()
              if v.get("paired") and v.get("pdbs", {}).get("unphos")
              and v.get("pdbs", {}).get("phos")]
    PURE = set()
GENE = {k: v["gene"].lower() for k, v in BENCHMARK_PROTEINS.items()}
SHORT = {
    "CDK2":"CDK2", "ERK2":"ERK2", "Aurora_A":"AURKA", "IRF3":"IRF3",
    "P38alpha":"p38", "CDK7":"CDK7", "c-Src_kinase":"Src", "AMPK_alpha2":"AMPK",
    "Glycogen_Phosphorylase":"PYGM", "CDK1":"CDK1", "ERK1":"ERK1",
    "GSK3B":"GSK3B", "PAK1":"PAK1", "AKT1":"AKT1", "TBK1":"TBK1", "IKBKB":"IKBb",
}
for _k, _v in BENCHMARK_PROTEINS.items():
    SHORT.setdefault(_k, _v["gene"])

def _sites(key):
    """从 benchmark 自动读取目标磷酸化位点（UniProt 坐标）"""
    v = BENCHMARK_PROTEINS.get(key, {})
    return [int("".join(c for c in s if c.isdigit())) for s in v.get("phos_sites", [])]
SITE = {k: _sites(k) for k in BENCHMARK_PROTEINS}

# B3：记录链回退事件，供全蛋白一致性审计
CHAIN_LOG = []


def ca_from_cif(path):
    """返回 (coords[n,3], seq[str], resids[list])；纳入 SEP/TPO/PTR（B2）。"""
    coords, seq, resids, meta = load_ca(path, chain=None)
    CHAIN_LOG.append({"file": Path(path).name, "side": "af3", **meta})
    return coords, seq, resids


def ca_from_pdb(path, chain="A"):
    """返回 (coords[n,3], seq, resids)；纳入 SEP/TPO/PTR（B2）；
    链缺失时回退并记录（B3）。"""
    coords, seq, resids, meta = load_ca(path, chain=chain)
    CHAIN_LOG.append({"file": Path(path).name, "side": "pdb", "requested": chain, **meta})
    return coords, seq, resids


def align_pairs(seqA, seqB):
    """序列对齐，返回 [(idxA, idxB)] 的匹配对（仅高置信匹配块）。"""
    sm = SequenceMatcher(None, seqA, seqB, autojunk=False)
    pairs = []
    for blk in sm.get_matching_blocks():
        for k in range(blk.size):
            pairs.append((blk.a + k, blk.b + k))
    return pairs


def kabsch(P, Q):
    """求把 P 叠合到 Q 的 (R, t)：P@R + t ≈ Q。

    注意：本机 numpy 以 Apple Accelerate 作 BLAS，其 (3×N)·(N×3) 核会
    对**有限输入**伪设 divide/overflow/underflow FPU 标志，numpy 随即在
    matmul 处报告 RuntimeWarning（已核实数值与纯 numpy einsum 完全一致）。
    故此处统一在 errstate 中屏蔽这四个类别；不影响数值。
    """
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        H = Pc.T @ Qc                   # 3×3 协方差（Accelerate 此处可能伪报警）
        U, _, Vt = np.linalg.svd(H)
        det = float(np.linalg.det(U @ Vt))
        d = 1.0 if det >= 0 else -1.0
        D = np.diag([1.0, 1.0, d])
        R = U @ D @ Vt                  # 行向量约定：P @ R ≈ Q
        t = Q.mean(0) - P.mean(0) @ R
    return R, t


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def af3_state_files(gene, state):
    """列出某态全部可用 AF3 模型文件（3 seeds × 5 models = 15）。

    List all available AF3 model files for one state (3 seeds × 5 models = 15).
    """
    d = CIF / gene / state
    if not d.exists():
        return []
    return sorted(d.glob("fold_*_seed*_model_*.cif"))


def compare(protein, model_idx=None):
    """以 PDB unphos 为共同框架，用序列对齐建立残基对应后比较位移。

    model_idx=None → 聚合全部可用模型（seed×model，通常 15 个）并取均值；
    model_idx=k    → 仅用第 k 号 model（跨所有 seed）。
    每蛋白实际独立模型数在结果字段 `n_models` 中显式给出（B4）。

    With model_idx=None all available models are aggregated (mean over models);
    with model_idx=k only model index k is used. The exact per-protein model
    count is reported in `n_models` (B4).
    """
    gene = GENE[protein]
    info = BENCHMARK_PROTEINS[protein]
    pdbs = info["pdbs"]
    chain = info.get("chain", "A")

    # ---- AF3 模型清单（B4：不再只取单个 CIF，而是遍历 seed×model）----
    ph_files = af3_state_files(gene, "phos")
    up_files = af3_state_files(gene, "unphos")
    if model_idx is not None:
        tag = f"_model_{model_idx}.cif"
        ph_files = [f for f in ph_files if f.name.endswith(tag)]
        up_files = [f for f in up_files if f.name.endswith(tag)]
    n_pairs = min(len(ph_files), len(up_files))
    if n_pairs == 0:
        return None

    # ---- PDB 侧（与模型无关，仅计算一次）----
    pdb_chains = info.get("pdb_chains", {})
    up_id = next((x for x in pdbs.get("unphos", []) if (PDBDIR/f"{x.upper()}.pdb").exists()), None)
    ph_id = next((x for x in pdbs.get("phos", []) if (PDBDIR/f"{x.upper()}.pdb").exists()), None)
    if not (up_id and ph_id):
        return None
    # 每个结构单独取链（CDK7: 7B5O 在链 J、8PYR 在链 A）
    ch_up = pdb_chains.get(up_id, pdb_chains.get(up_id.upper(), chain))
    ch_ph = pdb_chains.get(ph_id, pdb_chains.get(ph_id.upper(), chain))
    r_unph = ca_from_pdb(PDBDIR/f"{up_id.upper()}.pdb", ch_up)
    r_phos = ca_from_pdb(PDBDIR/f"{ph_id.upper()}.pdb", ch_ph)
    same_structure = (up_id.upper() == ph_id.upper())

    sites = SITE.get(protein, [])
    per_model, first_per_res = [], []

    for k in range(n_pairs):
        a_phos = ca_from_cif(ph_files[k])
        a_unph = ca_from_cif(up_files[k])

        # ---- 四方序列对齐（AF3 phos/unphos 与 PDB phos/unphos 长度可能不同） ----
        base_seq = a_unph[1]                                   # 以 AF3 unphos 序列为基准
        m_aph = dict(align_pairs(base_seq, a_phos[1]))         # af3_unphos idx -> af3_phos idx
        m_up  = dict(align_pairs(base_seq, r_unph[1]))         # -> pdb_unphos idx
        m_ph  = dict(align_pairs(base_seq, r_phos[1])) if not same_structure else m_up
        common = sorted(set(m_aph) & set(m_up) & set(m_ph))
        if len(common) < 30:
            continue

        ia_un = np.array(common)
        ia_ph = np.array([m_aph[i] for i in common])
        ip_up = np.array([m_up[i] for i in common])
        ip_ph = np.array([m_ph[i] for i in common])
        resids = [a_unph[2][i] for i in common]

        # ---- 叠合到 PDB unphos 框架 ----
        # 整个变换块置于 errstate 中：Accelerate BLAS 的 (N×3)·(3×3) 会伪报
        # FPU 标志（同 kabsch 注释），此处一并屏蔽，不影响数值。
        with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
            Q = r_unph[0][ip_up]
            if same_structure:
                pdb_phos_al = r_phos[0][ip_ph]
            else:
                P = r_phos[0][ip_ph]
                R0, t0 = kabsch(P, Q)
                pdb_phos_al = P @ R0 + t0

            # 关键：AF3 每次预测使用任意坐标系，两态必须各自独立叠合到同一参考框架，
            # 否则 d_pred 会混入整体刚体位移（曾导致 25-39 Å 的假信号）。
            P_au = a_unph[0][ia_un]
            R_u, t_u = kabsch(P_au, Q)
            af3_unphos_al = P_au @ R_u + t_u
            P_ap = a_phos[0][ia_ph]
            R_p, t_p = kabsch(P_ap, Q)
            af3_phos_al = P_ap @ R_p + t_p

            d_exp = pdb_phos_al - Q
            d_pred = af3_phos_al - af3_unphos_al
            mag_exp = np.linalg.norm(d_exp, axis=1)
            mag_pred = np.linalg.norm(d_pred, axis=1)

            # 投影约定敏感性（W1）：直接约定 = 两个 AF3 态先互相叠合再取差
            # （度量模型自身构象变化，不参照实验坐标；幅值与坐标系无关）。
            R_d, t_d = kabsch(P_ap, P_au)
            d_pred_direct = (P_ap @ R_d + t_d) - P_au
            mag_pred_direct = np.linalg.norm(d_pred_direct, axis=1)

        res = {
            "n_common": len(common),
            "rmsd_af3unphos_vs_pdb": round(rmsd(af3_unphos_al, Q), 2),
            "mean_disp_exp": round(float(mag_exp.mean()), 2),
            "mean_disp_pred": round(float(mag_pred.mean()), 2),
            "mean_disp_pred_direct": round(float(mag_pred_direct.mean()), 2),
            "max_disp_exp": round(float(mag_exp.max()), 2),
            "max_disp_pred": round(float(mag_pred.max()), 2),
        }

        if not same_structure:
            if mag_exp.std() > 1e-6 and mag_pred.std() > 1e-6:
                res["pearson_magnitude"] = round(float(np.corrcoef(mag_exp, mag_pred)[0, 1]), 3)
            sel = mag_exp > 1.0
            res["mean_cosine"] = None          # 默认无值，避免下游 KeyError
            res["n_dir_eval"] = int(sel.sum())
            if sel.sum() >= 5:
                u = d_exp[sel] / np.linalg.norm(d_exp[sel], axis=1, keepdims=True)
                vn = np.linalg.norm(d_pred[sel], axis=1, keepdims=True)
                v = d_pred[sel] / np.maximum(vn, 1e-9)
                cos = (u * v).sum(1)
                res["mean_cosine"] = round(float(cos.mean()), 3)
                res["frac_cosine_pos"] = round(float((cos > 0).mean()), 3)
            k20 = min(20, len(common))
            top_e = set(np.array(resids)[np.argsort(-mag_exp)[:k20]])
            top_p = set(np.array(resids)[np.argsort(-mag_pred)[:k20]])
            res["top20_overlap"] = len(top_e & top_p)

        if sites:
            near = [i for i, r in enumerate(resids) if any(abs(r - s) <= SITE_WIN for s in sites)]
            if len(near) >= 3:
                rest = [i for i in range(len(resids)) if i not in set(near)]
                res["site_window_n"] = len(near)
                res["site_mean_disp_exp"] = round(float(mag_exp[near].mean()), 2)
                res["site_mean_disp_pred"] = round(float(mag_pred[near].mean()), 2)
                if rest:
                    res["rest_mean_disp_exp"] = round(float(mag_exp[rest].mean()), 2)
                    res["rest_mean_disp_pred"] = round(float(mag_pred[rest].mean()), 2)

        per_model.append(res)
        if not first_per_res:      # 逐残基谱仅取首个模型，避免 15× 冗余
            first_per_res = [{"protein": protein, "short": SHORT[protein], "resid": int(r),
                              "disp_exp": round(float(mag_exp[i]), 3),
                              "disp_pred": round(float(mag_pred[i]), 3),
                              "is_site_window": any(abs(r - s) <= SITE_WIN for s in sites)}
                             for i, r in enumerate(resids)]

    if not per_model:
        return {"protein": protein, "short": SHORT[protein], "n_common": 0,
                "skipped": "insufficient common residues"}, []

    # ---- 跨模型聚合（B4）----
    agg = {"protein": protein, "short": SHORT[protein],
           "unphos_pdb": up_id, "phos_pdb": ph_id, "same_structure": same_structure,
           "unphos_chain": ch_up, "phos_chain": ch_ph,
           "pure_pair": (protein in PURE), "n_models": len(per_model)}
    INT_KEYS = ("n_common", "n_dir_eval", "site_window_n", "top20_overlap")
    num_keys = sorted({k for m in per_model for k, v in m.items()
                       if isinstance(v, (int, float)) and not isinstance(v, bool) and v is not None})
    for key in num_keys:
        vals = [m[key] for m in per_model if m.get(key) is not None]
        if not vals:
            continue
        mv = float(np.mean(vals))
        agg[key] = int(round(mv)) if key in INT_KEYS else round(mv, 3)
    # 模型间离散（用于 Fig2 误差棒；种子/模型不是生物学重复）
    # 同时报告模型内不稳定度（W10.iv）：同一蛋白 15 个模型的预测位移
    # min/max/CV，用于说明"近常量预测位移"在各蛋白内的稳健程度。
    dvals = [m["mean_disp_pred"] for m in per_model if m.get("mean_disp_pred") is not None]
    if len(dvals) > 1:
        dmean = float(np.mean(dvals))
        agg["sd_mean_disp_pred"] = round(float(np.std(dvals, ddof=1)), 3)
        agg["model_disp_min"] = round(float(np.min(dvals)), 3)
        agg["model_disp_max"] = round(float(np.max(dvals)), 3)
        agg["model_disp_cv"] = (round(float(np.std(dvals, ddof=1) / dmean), 3)
                                if dmean > 1e-9 else None)
    return agg, first_per_res


def bootstrap_ci(x, y, n_boot=10000, seed=42):
    """Pearson r 的 bootstrap 95% CI；重采样单位 = 蛋白（pair），非模型。

    Bootstrap 95% CI for Pearson r; the resampling unit is the **protein**
    (a paired observation), not the model, since seeds are not biological
    replicates (B4).
    """
    n = len(x)
    if n < 4 or np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return None, None
    rng = np.random.default_rng(seed)
    rs = []
    for _ in range(n_boot):
        s = rng.integers(0, n, n)
        xs, ys = x[s], y[s]
        if np.std(xs) < 1e-12 or np.std(ys) < 1e-12:
            continue
        rs.append(float(np.corrcoef(xs, ys)[0, 1]))
    if not rs:
        return None, None
    return float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))


def cohort_stats(results):
    """队列级核心统计：预测位移 vs 实验位移的 Pearson r + bootstrap CI（B4）。

    Cohort-level core statistic: Pearson r between per-protein mean predicted
    and experimental displacement, with bootstrap CI and leave-one-out.
    仅纳入有真实实验位移的蛋白（两态为不同 PDB 结构）。
    """
    sub = [r for r in results
           if not r.get("same_structure")
           and r.get("mean_disp_exp") is not None
           and r.get("mean_disp_pred") is not None]
    if len(sub) < 4:
        return None
    exp = np.array([r["mean_disp_exp"] for r in sub], float)
    pred = np.array([r["mean_disp_pred"] for r in sub], float)

    r = float(np.corrcoef(exp, pred)[0, 1])
    lo, hi = bootstrap_ci(exp, pred)
    slope = float(np.polyfit(exp, pred, 1)[0]) if np.std(exp) > 1e-12 else None

    try:
        from scipy import stats as _st
        rho, sp = _st.spearmanr(exp, pred)
        rho, sp = float(rho), float(sp)
    except Exception:
        rho, sp = None, None

    # leave-one-out：识别杠杆点
    loo = []
    for i in range(len(sub)):
        m = np.ones(len(sub), bool); m[i] = False
        if np.std(exp[m]) < 1e-12 or np.std(pred[m]) < 1e-12:
            continue
        loo.append((sub[i]["short"], float(np.corrcoef(exp[m], pred[m])[0, 1])))

    out = {
        "n_proteins": len(sub),
        "pearson_r": round(r, 4),
        "pearson_CI95": [None if lo is None else round(lo, 4),
                         None if hi is None else round(hi, 4)],
        "slope_pred_per_exp": None if slope is None else round(slope, 4),
        "spearman_rho": None if rho is None else round(rho, 4),
        "spearman_p": None if sp is None else round(sp, 4),
        "mean_pred": round(float(pred.mean()), 3),
        "sd_pred": round(float(pred.std(ddof=1)), 3),
        "range_exp": [round(float(exp.min()), 2), round(float(exp.max()), 2)],
        "leave_one_out": [{"removed": k, "r": round(v, 4)} for k, v in loo],
    }
    if loo:
        worst = max(loo, key=lambda kv: abs(kv[1] - r))
        out["max_leverage"] = {"removed": worst[0], "r_without": round(worst[1], 4)}

    # P3 覆盖失败对敏感性（W3）：剔除两个 P3 失败对（AMPK_alpha2、TBK1）
    # 这两个配对未通过 P3 覆盖门（<0.5），但仍被计入主队列，故给出剔除后的稳健性。
    # Sensitivity excluding the two P3-coverage failures (AMPK_alpha2, TBK1),
    # which were admitted to the primary cohort despite failing the P3 gate.
    EXCL = {"AMPK_alpha2", "TBK1"}
    sub2 = [r for r in sub if r.get("protein") not in EXCL]
    if len(sub2) >= 4:
        e2 = np.array([r["mean_disp_exp"] for r in sub2], float)
        p2 = np.array([r["mean_disp_pred"] for r in sub2], float)
        if np.std(e2) > 1e-12 and np.std(p2) > 1e-12:
            r2 = float(np.corrcoef(e2, p2)[0, 1])
            lo2, hi2 = bootstrap_ci(e2, p2)
            out["excl_p3_failures"] = {
                "removed": sorted(EXCL), "n_proteins": len(sub2),
                "pearson_r": round(r2, 4),
                "pearson_CI95": [None if lo2 is None else round(lo2, 4),
                                 None if hi2 is None else round(hi2, 4)],
            }

    # 纯度分层不确定性（R5）：pure / impure 两层各自的 r 与 bootstrap 95% CI，
    # 并用置换检验给出两层差异的 p 值（层标签随机重排 20k 次）。
    # Purity-stratified uncertainty (R5): each stratum's r with a bootstrap 95% CI,
    # plus a permutation p-value for the difference (20k label reshuffles).
    strata = {}
    for name, want in (("pure", True), ("impure", False)):
        s2 = [r for r in sub if bool(r.get("pure_pair")) == want]
        if len(s2) < 4:
            continue
        e2 = np.array([r["mean_disp_exp"] for r in s2], float)
        p2 = np.array([r["mean_disp_pred"] for r in s2], float)
        if np.std(e2) < 1e-12 or np.std(p2) < 1e-12:
            continue
        l2, h2 = bootstrap_ci(e2, p2)
        strata[name] = {
            "n_proteins": len(s2),
            "pearson_r": round(float(np.corrcoef(e2, p2)[0, 1]), 4),
            "pearson_CI95": [None if l2 is None else round(l2, 4),
                             None if h2 is None else round(h2, 4)],
        }
    if set(strata) == {"pure", "impure"}:
        flag = np.array([bool(r.get("pure_pair")) for r in sub])
        obs = abs(strata["pure"]["pearson_r"] - strata["impure"]["pearson_r"])
        rng = np.random.default_rng(7)
        hits = 0
        for _ in range(20000):
            perm = rng.permutation(flag)
            a = perm if perm.dtype == bool else perm.astype(bool)
            ea, pa = exp[~a], pred[~a]
            eb, pb = exp[a], pred[a]
            if min(len(ea), len(eb)) < 4:
                continue
            if min(np.std(ea), np.std(pa), np.std(eb), np.std(pb)) < 1e-12:
                continue
            d = abs(np.corrcoef(ea, pa)[0, 1] - np.corrcoef(eb, pb)[0, 1])
            if d >= obs:
                hits += 1
        strata["difference_permutation_p"] = round(hits / 20000, 4)
    if strata:
        out["purity_strata"] = strata
    return out


def projection_sensitivity(results):
    """两种投影约定下的队列级统计（W1 敏感性分析）。

    Cohort-level statistics under two projection conventions:
      indep  — both AF3 states independently superposed onto the PDB reference
               (the protocol the manuscript reports);
      direct — the two AF3 states superposed onto each other (the model's own
               internal phospho-induced change, no experimental reference).
    A "shared-rigid-transform" third convention is not defined here because the
    two AF3 states are independent predictions in unrelated coordinate frames.
    All use shipped models only; no new predictions.
    """
    def _one(predkey):
        sub = [r for r in results
               if not r.get("same_structure")
               and r.get("mean_disp_exp") is not None
               and r.get(predkey) is not None]
        if len(sub) < 4:
            return {"n_proteins": len(sub), "note": "insufficient proteins"}
        exp = np.array([r["mean_disp_exp"] for r in sub], float)
        pred = np.array([r[predkey] for r in sub], float)
        r = float(np.corrcoef(exp, pred)[0, 1])
        lo, hi = bootstrap_ci(exp, pred)
        return {
            "n_proteins": len(sub),
            "pearson_r": round(r, 4),
            "pearson_CI95": [None if lo is None else round(lo, 4),
                             None if hi is None else round(hi, 4)],
            "mean_pred": round(float(pred.mean()), 3),
            "sd_pred": round(float(pred.std(ddof=1)), 3),
        }

    return {
        "description": ("同一批已生成 AF3 模型在两种叠合/投影约定下的队列级结果；"
                        "用于证明结论不依赖于叠合约定（W1）。"),
        "protocols": {
            "indep": _one("mean_disp_pred"),
            "direct": _one("mean_disp_pred_direct"),
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", type=int, default=None,
                    help="仅用第 k 号 model（跨所有 seed）；默认聚合每蛋白全部模型")
    args = ap.parse_args()

    results, all_per = [], []
    print("=" * 92)
    print("  AF3 预测位移 vs PDB 实验位移（B4：默认聚合每蛋白全部 seed×model）")
    print("=" * 92)
    for p in PAIRED:
        try:
            out = compare(p, args.model)
        except Exception as e:
            print(f"  {p:<26} 失败: {type(e).__name__}: {e}")
            continue
        if not out:
            print(f"  {p:<26} 数据不足，跳过")
            continue
        r, per = out
        if r.get("skipped"):
            print(f"  {r['short']:<8} 跳过: {r['skipped']} (n={r['n_common']})")
            continue
        results.append(r); all_per.extend(per)
        nm = r.get("n_models", 1)
        if r["same_structure"]:
            print(f"  {r['short']:<8} n={r['n_common']:<4} models={nm:<3} "
                  f"同一PDB({r['unphos_pdb']})，无实验位移 | AF3位移均值 {r['mean_disp_pred']:.2f} Å")
        else:
            rho = r.get("pearson_magnitude")
            cos = r.get("mean_cosine")
            rho_s = f"{rho:+.3f}" if rho is not None else "  n/a"
            cos_s = f"{cos:+.3f}" if cos is not None else "  n/a"
            print(f"  {r['short']:<8} n={r['n_common']:<4} models={nm:<3} "
                  f"{r['unphos_pdb']}/{r['phos_pdb']:<5} r={rho_s} cos={cos_s} "
                  f"top20∩={r.get('top20_overlap')}/20 | 位移 exp {r['mean_disp_exp']:.2f} "
                  f"vs pred {r['mean_disp_pred']:.2f} Å")

    stats = cohort_stats(results)
    if stats:
        ci = stats["pearson_CI95"]
        ci_s = (f"[{ci[0]:+.4f}, {ci[1]:+.4f}]" if ci[0] is not None and ci[1] is not None else "n/a")
        rho_s = (f"{stats['spearman_rho']:+.4f} (p={stats['spearman_p']:.3f})"
                 if stats["spearman_rho"] is not None else "n/a")
        print("\n" + "-" * 92)
        print(f"  队列核心统计（以蛋白为独立单位，n={stats['n_proteins']}；模型数非生物学重复）")
        print(f"    Pearson r = {stats['pearson_r']:+.4f}   95% CI {ci_s}")
        print(f"    slope = {stats['slope_pred_per_exp']:+.4f} Å/Å   Spearman ρ = {rho_s}")
        print(f"    AF3 预测位移 {stats['mean_pred']:.2f} ± {stats['sd_pred']:.2f} Å；"
              f"实验位移跨度 {stats['range_exp'][0]:.2f}–{stats['range_exp'][1]:.2f} Å")
        if stats.get("max_leverage"):
            ml = stats["max_leverage"]
            print(f"    最大杠杆点：剔除 {ml['removed']} 后 r = {ml['r_without']:+.4f}")
        print("-" * 92)

    (OUT/"af3_vs_pdb_summary.json").write_text(
        json.dumps({"per_protein": results,
                    "cohort_stats": stats,
                    "model_mode": "all" if args.model is None else f"model_{args.model}"},
                   indent=2, ensure_ascii=False))
    proj = projection_sensitivity(results)
    (OUT/"projection_sensitivity.json").write_text(
        json.dumps(proj, indent=2, ensure_ascii=False))
    print("\n  投影约定敏感性（W1）:")
    for name, s in proj["protocols"].items():
        if "pearson_r" in s:
            print(f"    {name:<7} n={s['n_proteins']}  mean_pred={s['mean_pred']:.3f} Å  "
                  f"r={s['pearson_r']:+.4f}  CI={s['pearson_CI95']}")
    print(f"[OK] {OUT}/projection_sensitivity.json")
    if all_per:
        import csv
        with open(OUT/"af3_vs_pdb_per_residue.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_per[0].keys())); w.writeheader(); w.writerows(all_per)
    if CHAIN_LOG:      # B3：链回退审计
        (OUT/"chain_resolution_log_compare.json").write_text(
            json.dumps(CHAIN_LOG, indent=2, ensure_ascii=False))
    print(f"\n[OK] {OUT}/af3_vs_pdb_summary.json")
    return results


if __name__ == "__main__":
    main()
