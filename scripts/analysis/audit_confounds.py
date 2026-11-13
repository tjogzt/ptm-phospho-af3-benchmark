#!/usr/bin/env python3
"""
audit_confounds.py — 配对结构共变量审计（W4）
================================================
Confound audit for the phospho/unphospho pairs (W4).

中文（Chinese）
---------------
P5 目前只检查「链组成是否一致」，并不检查两态之间是否还存在**其它**共变量。
本脚本直接从 PDB 记录（TITLE / HET / HETATM / MODRES / SEQADV / DBREF）抽取这些
共变量，逐蛋白落盘，用于把正文的 "pure" 收紧为可辩护的 "chain-composition-matched"：

  - ligand                两态结合的配体（抑制剂 / ATP / 类似物）名称
  - mutation              点突变（来自 SEQADV，如 AURKA 4J8M 的 T287D、PAK1 1YHW 的 K299R）
  - phos_modres_sites_pdb  目标蛋白链上、MODRES 记录的磷酸位点（**PDB 编号**）
  - phos_modres_sites_uniprot 同一批位点经序列比对映射到 **UniProt 编号**
  - modelled_sites        AF3 建模时实际指定的位点（UniProt 编号，来自冻结 manifest）
  - extra_phos_sites      实验 phos 态存在、但 AF3 未建模的额外磷酸位点（UniProt 编号）
  - af3_input             AF3 的输入状态（恒为 apo 全长序列，不含任何配体）

编号约定（R3 修复 / numbering convention, R3 fix）
--------------------------------------------------
早期版本直接把 MODRES 的 PDB 残基号与 manifest 的 UniProt 残基号相减，导致两类
伪差异：(i) 物种/成熟链造成的整体编号偏移（ERK2 2ERK 为鼠源 P63086，Thr183/Tyr185
对应人源 Thr185/Tyr187；PYGM 1LWO 以成熟蛋白编号，Ser14 = UniProt Ser15）；
(ii) 不同亚基链（AMPK 6B1U 的 SEP B 108 位于 β1 亚基 Q9Y478，并非 α2 的位点）。
本版本对每条链做「链序列 → 目标 UniProt 序列」的序列比对；仅把**与目标蛋白同源**
（比对覆盖 ≥60% 的链）上的 MODRES 磷酸位点计入，并映射到 UniProt 编号后再比较。

注意：AF3 在**apo 全长序列**上运行，因此晶体配体/突变诱导的构象差异在预测
侧根本不可见；这正是"配对纯度"必须在**实验配对**层面审计的原因。

English
-------
P5 checks chain composition only. This script extracts the *other* covariates
(ligands, point mutations, the full MODRES set) directly from the shipped PDB
headers so the paper can narrow "pure" to the defensible
"chain-composition-matched". MODRES site numbers are mapped from PDB to UniProt
coordinates by sequence alignment, and only chains homologous to the target
protein are considered, so numbering offsets and hetero-subunit phosphosites are
not misreported as extra sites.

输出 / Output：
  results/analysis/pair_confound_audit.json
用法 / Usage：
  python3 scripts/analysis/audit_confounds.py
"""
import json
import sys
from difflib import SequenceMatcher
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from struct_utils import MODIFIED_AA          # 氨基酸类修饰残基白名单

PDBDIR = _ROOT / "data" / "pdb"
OUT = _ROOT / "results" / "analysis"
MANIFEST = _ROOT / "results" / "manifest_frozen.json"
UNIPROT_JSON = _ROOT / "data" / "curated" / "uniprot_sequences.json"

PHOS_RES = ("SEP", "TPO", "PTR")

# 常见溶剂/离子/缓冲成分，不计作生物学配体
IGNORE_HET = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "SO4", "PO4", "ACT", "ACY", "FMT",
    "MES", "TRS", "EPE", "IMD", "NA", "CL", "K", "MG", "CA", "ZN", "MN",
    "FE", "CU", "NI", "CD", "BR", "IOD", "NO3", "CIT", "TLA", "MPD", "DMS",
    "BME", "DTT", "GSH", "FLC", "PGE", "1PE", "P6G", "XPE", "BCT", "CO3",
    "CAF",  # 4J8M 的 S-二甲基砷酰半胱氨酸，非生物学配体
}

_THREE2ONE = {  # 含修饰残基的三字母 → 单字母
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C",
    "MLY": "K", "KCX": "K", "HYP": "P", "PCA": "Q",
}


def parse_pdb_meta(path):
    """从 PDB 头/记录抽取配体、修饰残基（按链）与点突变。"""
    title_parts, het, modres, modres_chain, seqadv = [], {}, {}, {}, []
    dbref = {}
    if not path.exists():
        return None
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "TITLE":
                title_parts.append(line[10:].rstrip())
            elif rec == "HET":
                name = line[7:10].strip().upper()
                if name and name not in IGNORE_HET:
                    het[name] = het.get(name, 0) + 1
            elif rec == "HETATM":
                name = line[17:20].strip().upper()
                if name and name not in IGNORE_HET:
                    het.setdefault(name, 0)
            elif rec == "MODRES":
                rn = line[12:15].strip().upper()
                ch = line[16]
                num = line[18:22].strip()
                if rn:
                    modres.setdefault(rn, []).append(num)
                    modres_chain.setdefault(ch, {})[num] = rn
            elif rec == "DBREF":
                ch = line[12]
                dbref[ch] = line[33:41].strip()
            elif rec == "SEQADV":
                seqadv.append(line[10:].strip())
    ligands = sorted(n for n in het if n not in MODIFIED_AA)
    mod_aa = sorted(n for n in het if n in MODIFIED_AA)
    return {
        "title": " ".join(title_parts).strip(),
        "ligands": ligands,
        "modified_aa": mod_aa,
        "modres": {k: sorted(set(v)) for k, v in sorted(modres.items())},
        "modres_chain": {c: dict(sorted(d.items())) for c, d in sorted(modres_chain.items())},
        "dbref": dbref,
        "seqadv": seqadv,
    }


def chain_ca_seq(path):
    """各链的 Cα 残基序列 → {chain: [(resnum, one_letter), ...]}（按出现顺序，去重）。

    仅取 ATOM/HETATM 的 CA 原子，避免对同一残基重复计数。
    """
    chains = {}
    seen = set()
    if not path.exists():
        return chains
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            ch = line[21]
            key = (ch, line[22:27])
            if key in seen:
                continue
            aa = _THREE2ONE.get(line[17:20].strip().upper())
            if aa is None:
                continue
            try:
                num = int(line[22:26])
            except ValueError:
                continue
            seen.add(key)
            chains.setdefault(ch, []).append((num, aa))
    return chains


def map_pdb_to_uniprot(pairs, uniprot_seq):
    """用序列比对得到 (映射, 匹配数)；映射为 {PDB 残基号 → UniProt 位置(1-based)}。

    以链的单字母序列对目标 UniProt 序列做比对，仅采用匹配块，
    因此可自动吸收物种差异与成熟链造成的整体/局部编号偏移。
    """
    if not pairs or not uniprot_seq:
        return {}, 0
    pdb_seq = "".join(a for _, a in pairs)
    sm = SequenceMatcher(None, pdb_seq, uniprot_seq, autojunk=False)
    mapping, matched = {}, 0
    for blk in sm.get_matching_blocks():
        for i in range(blk.size):
            if pdb_seq[blk.a + i] == uniprot_seq[blk.b + i]:
                mapping[pairs[blk.a + i][0]] = blk.b + i + 1   # UniProt 1-based
                matched += 1
    return mapping, matched


def mutations(meta, uniprot=None):
    """从 SEQADV 提取点突变（形如 'MUTATION' 或 'K->R' 的冲突注释）。

    仅保留注释在**目标蛋白**上的记录：同一晶体中的配对亚基（如 CDK1 的
    cyclin B1，P14635）也会带 SEQADV 记录，按 dbAccession 列过滤可避免把
    非目标链的工程突变误计为 has_mutation（W6 修复）。
    """
    out = []
    for s in meta.get("seqadv", []):
        if not ("MUTATION" in s.upper() or "->" in s or "VARIANT" in s.upper()):
            continue
        acc = s[19:28].strip()          # SEQADV 的 dbAccession 列
        if uniprot and acc and acc != uniprot:
            continue                    # 注释位于非目标链，不计入
        out.append(s)
    return out


def main():
    mf = json.loads(MANIFEST.read_text())
    upseq = json.loads(UNIPROT_JSON.read_text()) if UNIPROT_JSON.exists() else {}
    rows = []
    print(f"{'蛋白':<12}{'配体(up/ph)':>20}{'突变':>6}{'MODRES(PDB)':>18}"
          f"{'→UniProt':>12}{'建模位点':>12}{'额外(UniProt)':>14}")
    print("-" * 100)
    for p in mf["proteins"]:
        if not p.get("paired"):
            continue
        up_id = p.get("unphos_pdb")
        ph_id = p.get("phos_pdb")
        mu = parse_pdb_meta(PDBDIR / f"{str(up_id).upper()}.pdb") if up_id else None
        mp = parse_pdb_meta(PDBDIR / f"{str(ph_id).upper()}.pdb") if ph_id else None
        if not mu or not mp:
            continue

        target_seq = upseq.get(p.get("uniprot"), {}).get("sequence", "")
        chains = chain_ca_seq(PDBDIR / f"{str(ph_id).upper()}.pdb")

        # 找出与目标蛋白同源（比对覆盖 ≥60%）的链，逐链建立 PDB→UniProt 映射
        same_protein_chains, pdb2up, matched_frac = {}, {}, {}
        for ch, pairs in chains.items():
            if len(pairs) < 20:
                continue
            m, matched = map_pdb_to_uniprot(pairs, target_seq)
            cov = matched / len(pairs) if pairs else 0.0
            matched_frac[ch] = round(cov, 3)
            db = mp.get("dbref", {}).get(ch, "")
            if cov >= 0.6 or db == p.get("uniprot") or ch == p.get("phos_chain"):
                same_protein_chains[ch] = cov
                pdb2up.update(m)

        # 目标蛋白链上的 MODRES 磷酸位点（PDB 编号）
        modres_pdb = sorted({num for ch in same_protein_chains
                             for num, rn in mp.get("modres_chain", {}).get(ch, {}).items()
                             if rn in PHOS_RES})
        modres_up = sorted({pdb2up[int(n)] for n in modres_pdb
                            if n.isdigit() and int(n) in pdb2up})
        modelled = sorted({int(s) for s in (p.get("sites_uniprot") or [])})
        extra = sorted(set(modres_up) - set(modelled))

        # 其它亚基（非目标蛋白）上的磷酸位点，仅供信息，不计入 extra
        other_phos = {ch: sorted(num for num, rn in mp.get("modres_chain", {}).get(ch, {}).items()
                                 if rn in PHOS_RES)
                      for ch in mp.get("modres_chain", {})
                      if ch not in same_protein_chains
                      and any(rn in PHOS_RES for rn in mp["modres_chain"][ch].values())}

        lig_up, lig_ph = set(mu["ligands"]), set(mp["ligands"])
        mut_up, mut_ph = mutations(mu, p.get("uniprot")), mutations(mp, p.get("uniprot"))

        rec = {
            "key": p["key"],
            "gene": p.get("gene"),
            "uniprot": p.get("uniprot"),
            "unphos_pdb": up_id, "phos_pdb": ph_id,
            "phos_chain": p.get("phos_chain"),
            "pure_pair": p.get("pure_pair"),
            "in_primary_cohort": p.get("in_primary_cohort"),
            "ligands_unphos": sorted(lig_up),
            "ligands_phos": sorted(lig_ph),
            "ligand_differs": lig_up != lig_ph,
            "ligands_only_unphos": sorted(lig_up - lig_ph),
            "ligands_only_phos": sorted(lig_ph - lig_up),
            "mutation_unphos": mut_up,
            "mutation_phos": mut_ph,
            "has_mutation": bool(mut_up or mut_ph),
            "same_protein_chains": {c: matched_frac.get(c) for c in sorted(same_protein_chains)},
            "phos_modres_sites_pdb": modres_pdb,
            "phos_modres_sites_uniprot": modres_up,
            "af3_modelled_sites": modelled,
            "extra_phos_sites": extra,
            "extra_phos_sites_other_chains": other_phos,
            "dual_or_multi_phospho": len(modres_up) > 1,
            "af3_input": "apo full-length UniProt; no ligands modelled",
            "composition_matched": (
                p.get("pure_pair") and not lig_up.symmetric_difference(lig_ph)
                and not (mut_up or mut_ph)
            ),
        }
        rows.append(rec)

        lig_s = f"{','.join(sorted(lig_up)) or '-'} / {','.join(sorted(lig_ph)) or '-'}"
        mut_s = "Y" if (mut_up or mut_ph) else "-"
        print(f"{p['key']:<12}{lig_s[:20]:>20}{mut_s:>6}"
              f"{('-,'.join(modres_pdb) or '-'):>18}"
              f"{('-,'.join(map(str, modres_up)) or '-'):>12}"
              f"{('-,'.join(map(str, modelled)) or '-'):>12}"
              f"{('-,'.join(map(str, extra)) or '-'):>14}")

    n_pure = sum(1 for r in rows if r["pure_pair"])
    n_lig = sum(1 for r in rows if r["pure_pair"] and r["ligand_differs"])
    n_mut = sum(1 for r in rows if r["pure_pair"] and r["has_mutation"])
    n_comp = sum(1 for r in rows if r["composition_matched"])
    n_extra = sum(1 for r in rows if r["extra_phos_sites"])
    n_extra_other = sum(1 for r in rows if r["extra_phos_sites_other_chains"])
    summary = {
        "n_pairs_audited": len(rows),
        "n_pure_p5": n_pure,
        "n_pure_with_ligand_difference": n_lig,
        "n_pure_with_point_mutation": n_mut,
        "n_pairs_with_extra_phos": n_extra,
        "n_pairs_with_extra_phos_other_chains": n_extra_other,
        "n_composition_matched": n_comp,
        "note": ("'pure' (P5) 仅表示链组成一致；配体/点突变/额外磷酸位点是额外的、"
                 "P5 未覆盖的共变量。AF3 在 apo 全长序列上运行，故这些实验侧差异在"
                 "预测侧不可见。extra_phos_sites 已将 MODRES 位点由 PDB 编号映射到 "
                 "UniProt 编号，并限定于与目标蛋白同源的链（R3）。"
                 "composition_matched 为更严格的标记。"),
    }
    (OUT / "pair_confound_audit.json").write_text(
        json.dumps({"summary": summary, "per_protein": rows},
                   indent=2, ensure_ascii=False))
    print("-" * 100)
    print(f"  P5-pure {n_pure}；其中配体不同 {n_lig}，含点突变 {n_mut}，"
          f"含额外磷酸位点 {n_extra}（另有 {n_extra_other} 对含异源链磷酸位点）；"
          f"严格 composition-matched {n_comp}")
    print(f"[OK] {OUT}/pair_confound_audit.json")


if __name__ == "__main__":
    main()
