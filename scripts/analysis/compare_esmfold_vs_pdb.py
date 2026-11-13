#!/usr/bin/env python3
"""
compare_esmfold_vs_pdb.py — B7：ESMFold 预测位移 vs PDB 实验位移
======================================================================
第二条对比线。用与 compare_af3_vs_pdb.py **完全相同**的坐标系与度量协议，
把 ESMFold（phosphomimetic 策略）的预测位移与实验位移、以及 AF3 的预测位移
放在同一框架下比较，回答："AF3 的无响应是否为模型特异现象？"

方法 / Method（与 AF3 线一致）：
  1. 以 PDB unphos 为共同参考框架 Q
  2. PDB phos → 叠合到 Q → Δ_exp
  3. ESMFold 两态**各自独立**叠合到 Q（每次预测各有任意坐标系）：
     Δ_pred_esm = (ESMFold_phos fit→Q) − (ESMFold_unphos fit→Q)
     独立叠合消除整体刚体位移；仅用单一变换会把两态坐标系差异误当作磷酸化位移。
  4. 共同残基上：幅度 Pearson r、方向余弦、top-20 重合度、位点窗统计

ESMFold 每态仅 1 个确定性结构（无 seed×model），故不做跨模型聚合。

输出 / Output：
  results/analysis/esmfold_vs_pdb_summary.json
  results/analysis/esmfold_vs_pdb_per_residue.csv

用法 / Usage：python3 scripts/analysis/compare_esmfold_vs_pdb.py
"""
import sys, json, csv
from pathlib import Path
import numpy as np
from difflib import SequenceMatcher

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
sys.path.insert(0, str(_HERE))
from benchmark_dataset import BENCHMARK_PROTEINS
from struct_utils import load_ca

PDBDIR = _ROOT / "data" / "pdb"
ESMPDB = _ROOT / "results" / "esm" / "pdb"
OUT = _ROOT / "results" / "analysis"
OUT.mkdir(parents=True, exist_ok=True)

SITE_WIN = 10

MANIFEST = _ROOT / "results" / "manifest_frozen.json"
_mf = json.loads(MANIFEST.read_text())
PRIMARY = list(_mf["cohorts"]["primary"])
PURE = set(_mf["cohorts"]["pure"])

SHORT = {
    "CDK2": "CDK2", "ERK2": "ERK2", "Aurora_A": "AURKA", "IRF3": "IRF3",
    "P38alpha": "p38", "CDK7": "CDK7", "c-Src_kinase": "Src",
    "AMPK_alpha2": "AMPK", "Glycogen_Phosphorylase": "PYGM", "CDK1": "CDK1",
    "GSK3B": "GSK3B", "PAK1": "PAK1", "AKT1": "AKT1", "TBK1": "TBK1",
    "IKBKB": "IKBb",
}


def sites_of(info):
    return list(info.get("phos_sites_uniprot") or
                [int("".join(c for c in s if c.isdigit())) for s in info.get("phos_sites", [])])


def align_pairs(a, b):
    sm = SequenceMatcher(None, a, b, autojunk=False)
    out = []
    for blk in sm.get_matching_blocks():
        for k in range(blk.size):
            out.append((blk.a + k, blk.b + k))
    return out


def kabsch(P, Q):
    """叠合 P 到 Q：P@R + t ≈ Q。errstate 屏蔽 Apple Accelerate 伪 FPU 警告。"""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        H = Pc.T @ Qc
        U, _, Vt = np.linalg.svd(H)
        d = 1.0 if float(np.linalg.det(U @ Vt)) >= 0 else -1.0
        R = U @ np.diag([1.0, 1.0, d]) @ Vt
        t = Q.mean(0) - P.mean(0) @ R
    return R, t


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def load_esm(gene, state):
    """读 ESMFold 预测 PDB（缺失返回 None）。"""
    p = ESMPDB / f"{gene}_{state}_esmfold.pdb"
    if not p.exists():
        return None
    c, s, r, _ = load_ca(p, chain=None)
    return (c, s, r) if len(s) else None


def compare(protein):
    info = BENCHMARK_PROTEINS.get(protein)
    if not info or not info.get("paired"):
        return None
    gene = info["gene"].lower()
    esm_up = load_esm(gene, "unphos")
    esm_ph = load_esm(gene, "phos")
    if esm_up is None or esm_ph is None:
        return {"protein": protein, "short": SHORT.get(protein, gene), "skipped": "no ESMFold PDB"}

    pdbs = info["pdbs"]
    pdb_chains = info.get("pdb_chains", {})
    up_id = next((x for x in pdbs.get("unphos", []) if (PDBDIR / f"{x.upper()}.pdb").exists()), None)
    ph_id = next((x for x in pdbs.get("phos", []) if (PDBDIR / f"{x.upper()}.pdb").exists()), None)
    if not (up_id and ph_id):
        return {"protein": protein, "short": SHORT.get(protein, gene), "skipped": "missing PDB"}
    ch_up = pdb_chains.get(up_id, pdb_chains.get(up_id.upper(), info.get("chain", "A")))
    ch_ph = pdb_chains.get(ph_id, pdb_chains.get(ph_id.upper(), info.get("chain", "A")))
    r_up = load_ca(PDBDIR / f"{up_id.upper()}.pdb", chain=ch_up)[:3]
    r_ph = load_ca(PDBDIR / f"{ph_id.upper()}.pdb", chain=ch_ph)[:3]
    same_structure = (up_id.upper() == ph_id.upper())

    base_seq = esm_up[1]                       # ESMFold unphos 序列为基准（= AF3 序列）
    # ESMFold phos 由 unphos 在目标位点做 phosphomimetic 替换得到，长度不变、残基顺序
    # 一一对应，故用**按位置映射**而非序列比对：序列比对会在被替换的位点上少匹配 1 个
    # 残基，导致 ESM 线与 AF3 线的共同残基集不一致（进而使"实验位移"随线而变）。
    # ESMFold phos is built from unphos by positional phosphomimetic substitution, so a
    # positional 1:1 map (not sequence alignment) is the correct correspondence; this
    # keeps the common-residue set identical to the AF3 line.
    if len(esm_ph[1]) == len(base_seq):
        m_eph = {i: i for i in range(len(base_seq))}
    else:
        m_eph = dict(align_pairs(base_seq, esm_ph[1]))
    m_up = dict(align_pairs(base_seq, r_up[1]))
    m_ph = dict(align_pairs(base_seq, r_ph[1])) if not same_structure else m_up
    common = sorted(set(m_eph) & set(m_up) & set(m_ph))
    if len(common) < 30:
        return {"protein": protein, "short": SHORT.get(protein, gene),
                "skipped": f"insufficient common residues ({len(common)})"}

    ie_un = np.array(common)
    ie_ph = np.array([m_eph[i] for i in common])
    ip_up = np.array([m_up[i] for i in common])
    ip_ph = np.array([m_ph[i] for i in common])
    resids = [esm_up[2][i] for i in common]
    sites = sites_of(info)

    with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
        Q = r_up[0][ip_up]
        if same_structure:
            pdb_phos_al = r_ph[0][ip_ph]
        else:
            P = r_ph[0][ip_ph]
            R0, t0 = kabsch(P, Q)
            pdb_phos_al = P @ R0 + t0

        E_u = esm_up[0][ie_un]
        R_u, t_u = kabsch(E_u, Q)
        esm_up_al = E_u @ R_u + t_u
        E_p = esm_ph[0][ie_ph]
        R_p, t_p = kabsch(E_p, Q)
        esm_ph_al = E_p @ R_p + t_p

        d_exp = pdb_phos_al - Q
        d_pred = esm_ph_al - esm_up_al
        mag_exp = np.linalg.norm(d_exp, axis=1)
        mag_pred = np.linalg.norm(d_pred, axis=1)

    rec = {
        "protein": protein, "short": SHORT.get(protein, gene),
        "gene": gene, "unphos_pdb": up_id, "phos_pdb": ph_id,
        "same_structure": same_structure, "pure_pair": protein in PURE,
        "n_common": len(common),
        "rmsd_esmumphos_vs_pdb": round(rmsd(esm_up_al, Q), 2),
        "mean_disp_exp": round(float(mag_exp.mean()), 3),
        "mean_disp_pred": round(float(mag_pred.mean()), 3),
        "max_disp_exp": round(float(mag_exp.max()), 3),
        "max_disp_pred": round(float(mag_pred.max()), 3),
        "amplitude_ratio": round(float(mag_pred.mean() / mag_exp.mean()), 3)
        if mag_exp.mean() > 1e-9 else None,
    }
    if not same_structure:
        if mag_exp.std() > 1e-6 and mag_pred.std() > 1e-6:
            rec["pearson_magnitude"] = round(float(np.corrcoef(mag_exp, mag_pred)[0, 1]), 3)
        sel = mag_exp > 1.0
        rec["n_dir_eval"] = int(sel.sum())
        rec["mean_cosine"] = None
        if sel.sum() >= 5:
            u = d_exp[sel] / np.linalg.norm(d_exp[sel], axis=1, keepdims=True)
            vn = np.linalg.norm(d_pred[sel], axis=1, keepdims=True)
            cos = (u * (d_pred[sel] / np.maximum(vn, 1e-9))).sum(1)
            rec["mean_cosine"] = round(float(cos.mean()), 3)
            rec["frac_cosine_pos"] = round(float((cos > 0).mean()), 3)
        k20 = min(20, len(common))
        top_e = set(np.array(resids)[np.argsort(-mag_exp)[:k20]])
        top_p = set(np.array(resids)[np.argsort(-mag_pred)[:k20]])
        rec["top20_overlap"] = len(top_e & top_p)
    if sites:
        near = [i for i, r in enumerate(resids) if any(abs(r - s) <= SITE_WIN for s in sites)]
        if len(near) >= 3:
            rec["site_window_n"] = len(near)
            rec["site_mean_disp_exp"] = round(float(mag_exp[near].mean()), 3)
            rec["site_mean_disp_pred"] = round(float(mag_pred[near].mean()), 3)

    per_res = [{"protein": protein, "short": rec["short"], "resid": int(r),
                "disp_exp": round(float(mag_exp[i]), 3),
                "disp_pred": round(float(mag_pred[i]), 3),
                "is_site_window": bool(sites and any(abs(r - s) <= SITE_WIN for s in sites))}
               for i, r in enumerate(resids)]
    return rec, per_res


def bootstrap_ci(x, y, n_boot=10000, seed=42):
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
    return (float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))) if rs else (None, None)


def cohort_stats(rows):
    """队列级核心统计：预测位移 vs 实验位移的 Pearson r + bootstrap CI。

    与 compare_af3_vs_pdb.cohort_stats **完全同构**（相同的筛选、键名与
    leave-one-out 杠杆点识别），以便 make_tables.py 用同一条代码路径处理
    两条对比线（避免 ESM 键名不匹配导致 ESM 列静默丢失）。

    Identical in structure (filter, keys, leave-one-out leverage) to
    compare_af3_vs_pdb.cohort_stats, so make_tables.py can process both
    comparison lines through the same code path.
    """
    sub = [r for r in rows
           if not r.get("same_structure")
           and r.get("mean_disp_exp") is not None
           and r.get("mean_disp_pred") is not None]
    if len(sub) < 4:
        return {"n_proteins": len(sub), "note": "insufficient proteins"}
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

    # leave-one-out：识别杠杆点（剔除后 r 变化最大者）
    # Leave-one-out: identify the leverage point (largest change in r when removed)
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
    return out


def main():
    rows, per_res_all = [], []
    for key in PRIMARY:
        res = compare(key)
        if res is None:
            continue
        if isinstance(res, tuple):
            rec, pr = res
            rows.append(rec)
            per_res_all.extend(pr)
        else:
            rows.append(res)

    stats = cohort_stats(rows)
    summary = {
        "predictor": "ESMFold v1 (facebook/esmfold_v1, transformers EsmForProteinFolding)",
        "phospho_model": "phosphomimetic (pS/pT->Asp, pY->Glu); ESMFold cannot take CCD_ modifications",
        "device_note": "run locally on Apple Silicon (MPS) or CPU",
        "n_proteins": len([r for r in rows if r.get("mean_disp_exp") is not None]),
        "cohort_stats": stats,
        "per_protein": rows,
        "caveat": ("ESMFold 线使用 phosphomimetic 近似而非真实磷酸基团；此线用于判断"
                   "'无响应'是否为 AF3 特有，而非等价于真实磷酸化建模。"),
    }
    (OUT / "esmfold_vs_pdb_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    if per_res_all:
        with open(OUT / "esmfold_vs_pdb_per_residue.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["protein", "short", "resid", "disp_exp",
                                              "disp_pred", "is_site_window"])
            w.writeheader()
            w.writerows(per_res_all)

    print("=" * 92)
    print("  B7 ESMFold（phosphomimetic）vs PDB 实验位移")
    print("=" * 92)
    print(f"  {'蛋白':<8}{'共同残基':>8}{'实验位移':>10}{'预测位移':>10}"
          f"{'幅度比':>8}{'向量r':>8}{'方向cos':>9}{'top20':>7}")
    print("-" * 92)
    for r in sorted(rows, key=lambda x: (x.get("mean_disp_exp") is None, x.get("mean_disp_exp") or 0)):
        if r.get("mean_disp_exp") is None:
            print(f"  {r.get('short','?'):<8}  [skip] {r.get('skipped','')}")
            continue
        def f(k, d=8):
            v = r.get(k)
            return f"{v:>{d}.3f}" if isinstance(v, (int, float)) else f"{'n/a':>{d}}"
        print(f"  {r['short']:<8}{r['n_common']:>8}{f('mean_disp_exp',10)}"
              f"{f('mean_disp_pred',10)}{f('amplitude_ratio',8)}"
              f"{f('pearson_magnitude',8)}{f('mean_cosine',9)}"
              f"{r.get('top20_overlap','n/a'):>7}")
    print("-" * 92)
    print(f"  队列核心（n={stats.get('n_proteins')}）：")
    if "pearson_r" in stats:
        ci = stats["pearson_CI95"]
        ci_s = (f"[{ci[0]:+.4f}, {ci[1]:+.4f}]"
                if ci[0] is not None and ci[1] is not None else "n/a")
        print(f"    Pearson r = {stats['pearson_r']:+.4f}  95% CI {ci_s}")
        if stats.get("spearman_rho") is not None:
            print(f"    Spearman rho = {stats['spearman_rho']:+.4f} (p={stats['spearman_p']})")
        if stats.get("slope_pred_per_exp") is not None:
            print(f"    slope = {stats['slope_pred_per_exp']:+.4f} Å/Å")
        print(f"    预测位移 = {stats['mean_pred']} ± {stats['sd_pred']} Å；"
              f"实验位移跨度 {stats['range_exp']}")
        if "max_leverage" in stats:
            ml = stats["max_leverage"]
            print(f"    最大杠杆点 {ml['removed']}；剔除后 r = {ml['r_without']:+.4f}")
    print("-" * 92)
    print(f"[OK] {OUT}/esmfold_vs_pdb_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
