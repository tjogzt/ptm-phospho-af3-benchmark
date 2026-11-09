#!/usr/bin/env python3
"""
struct_utils.py — 结构 IO 共享工具（B2/B3 修复核心）
====================================================================
Shared structure-IO utilities (single source of truth for B2/B3).

中文（Chinese）
---------------
本模块解决两个已确认缺陷：
  (B2) Cα 提取使用 `res.id[0] == " "` 过滤，会把以 HETATM 记录出现的
       **修饰残基（SEP/TPO/PTR）** 一并丢弃——这不仅发生在 PDB 侧，
       也发生在 AF3 输出的 CIF 侧（AF3 把 TPO/SEP 写作 HETATM）。
       结果：磷酸化位点残基在「实验 vs 预测」两侧都被排除，度量失真。
  (B3) 不同脚本用不同的链选择/回退逻辑，导致同一蛋白出现矛盾数值。
       本模块统一：优先使用调用方指定的链；缺失时才回退到首链，
       并**显式返回是否发生回退**，供审计记录。

English
-------
This module fixes two confirmed defects:
  (B2) Cα extraction filtered on `res.id[0] == " "`, silently dropping
       modified residues (SEP/TPO/PTR) recorded as HETATM — on BOTH the
       experimental PDB side and the AF3 CIF side (AF3 writes TPO/SEP as
       HETATM). Consequently, the phosphosite residue was excluded on both
       sides of the experimental-vs-predicted comparison.
  (B3) Different scripts used different chain-selection / fallback logic,
       producing contradictory values for the same protein. This module
       unifies selection: use the requested chain; fall back to the first
       chain only if absent, and explicitly report whether a fallback
       occurred (for auditing).

用法 / Usage
------------
    from struct_utils import load_ca, resolve_chain, aa_from_resname
    coords, seq, resids, meta = load_ca("path.pdb", chain="A")
"""
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

# 20 标准氨基酸（three-letter → one-letter）
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# 需纳入链序列的氨基酸类修饰残基（含磷酸化残基）→ 母体氨基酸
# Amino-acid-like modified residues (incl. phospho-residues) → parent residue
MODIFIED_AA = {
    "SEP": "S", "TPO": "T", "PTR": "Y", "NEP": "H",   # 磷酸化 / phospho
    "MLY": "K", "M3L": "K", "KCX": "K",               # 甲基化/羧化 Lys
    "MSE": "M", "FME": "M",                          # 硒代/甲酰 Met
    "CSO": "C", "CME": "C", "CSD": "C", "OCS": "C",  # 氧化/半胱磺酸 Cys
    "HYP": "P", "PCA": "Q", "SEP_": "S",
    "CSX": "C", "CAS": "C", "PHD": "D", "MEN": "N",
    "TYS": "Y", "SME": "M", "LLP": "K",
}


def aa_from_resname(resname):
    """修饰残基 / 标准残基 → 单字母；非氨基酸类返回 None。
    Map (modified) residue name to one-letter code; None if not amino-acid-like."""
    r = resname.strip().upper()
    if r in AA3:
        return AA3[r]
    return MODIFIED_AA.get(r)


def _is_amino_res(res):
    """True 若残基为氨基酸（含修饰），无论 ATOM/HETATM。"""
    return res.id[0][0] in (" ", "H") and aa_from_resname(res.get_resname()) is not None


def resolve_chain(model, requested):
    """返回 (chain_id, fallback_used)。

    优先使用 requested；缺失时回退到首链。仅当**显式请求了某条链却未找到**
    时才将 fallback 置 True（`chain=None` 表示"任取首链"，属预期行为）。
    """
    ids = [c.get_id() for c in model]
    if not ids:
        return None, True
    if requested and requested in ids:
        return requested, False
    return ids[0], bool(requested)


_cif = MMCIFParser(QUIET=True)
_pdb = PDBParser(QUIET=True)


def load_ca(path, chain=None):
    """从 PDB 或 CIF 读 Cα 坐标。

    返回 / Returns:
      coords  : np.ndarray (n, 3)
      seq     : str（修饰残基映射回母体单字母）
      resids  : list[int]（res.id[1]）
      meta    : dict(chain_used, fallback, n_res, fmt)

    磷酸化等修饰残基（HETATM）**会被纳入**（B2）。
    """
    path = Path(path)
    parser = _cif if path.suffix.lower() in (".cif", ".mmcif") else _pdb
    s = parser.get_structure("x", str(path))
    model = next(iter(s))
    chain_used, fallback = resolve_chain(model, chain)
    coords, seq, resids = [], [], []
    if chain_used is not None:
        for res in model[chain_used]:
            if "CA" not in res:
                continue
            aa = aa_from_resname(res.get_resname()) if _is_amino_res(res) else None
            if aa is None:
                continue
            coords.append(res["CA"].coord.astype(float))
            seq.append(aa)
            resids.append(res.id[1])
    meta = {"chain_used": chain_used, "fallback": fallback, "n_res": len(resids),
            "fmt": path.suffix.lower()}
    arr = np.array(coords) if coords else np.zeros((0, 3))
    return arr, "".join(seq), resids, meta


def load_ca_first_chain(path):
    """读 CIF/PDB 中首个含残基的链（AF3 输出常用于此）。"""
    return load_ca(path, chain=None)
