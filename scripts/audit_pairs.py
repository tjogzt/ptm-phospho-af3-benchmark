#!/usr/bin/env python3
"""
audit_pairs.py — 基准配对质量筛查（P0/P1 前置门禁）
=====================================================
对 benchmark_dataset.py 中每个蛋白的每个 PDB 做四项硬核验，
防止 7ML6（非 CDK7）、1T2K（非磷酸化且异结构域）类错误再次入库。

核验项：
  C1 身份核验    PDB 头 TITLE/COMPND 的分子名需与声明基因/protein 名相符
  C2 真实磷酸化  phos 态结构必须含真实磷酸基团（MODRES SEP/TPO/PTR 或磷酸根配体）
  C3 序列一致    与 UniProt canonical 对齐，覆盖度与一致性达阈值
  C4 结构域重叠  phos 与 unphos 覆盖的残基范围须重叠且覆盖磷酸化位点

判定：PASS / WARN / FAIL，并给出可执行处置建议。

输出：
  audit/pair_audit.json     机器可读明细
  audit/pair_audit.md       人读报告
  stdout                    汇总表

用法：python3 scripts/audit_pairs.py [--min-cov 0.5] [--min-overlap 20]
"""
import sys, json, argparse, re
from pathlib import Path
from collections import OrderedDict

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

PDBDIR = _ROOT / "data" / "pdb"
OUT = _ROOT / "audit"; OUT.mkdir(exist_ok=True)

PHOS_RES = {"SEP", "TPO", "PTR", "NEP"}
PHOS_LIG = {"PO4", "SO4"}          # PO4 需进一步确认是否为磷酸根
AA3 = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
}


# ---------- 基础解析（直接读文本，比 BioPython 更可控） ----------
def read_header(path):
    """返回 (title, compnd_names[], keywords)"""
    title, compnd, mols = [], [], []
    try:
        with open(path, errors="ignore") as f:
            for line in f:
                rec = line[:6].strip()
                if rec == "TITLE":
                    title.append(line[10:].strip())
                elif rec == "COMPND":
                    compnd.append(line[10:].strip())
                elif rec in ("ATOM", "HETATM"):
                    break
                elif rec == "SEQRES" and not compnd:
                    break
    except Exception:
        pass
    t = " ".join(title)
    m = re.findall(r"MOLECULE:\s*([^;]+)", " ".join(compnd))
    return t, [x.strip() for x in m]


def read_modres(path):
    """返回 {chain: {resseq: (resname, std_res)}} 以及磷酸根配体计数"""
    mods, po4 = {}, 0
    try:
        with open(path, errors="ignore") as f:
            for line in f:
                rec = line[:6].strip()
                if rec == "MODRES":
                    resn = line[12:15].strip()
                    ch = line[16]
                    seq = line[18:22].strip()
                    std = line[24:27].strip()
                    if resn in PHOS_RES:
                        mods.setdefault(ch, {})[int(seq)] = (resn, std)
                elif rec == "HETATM":
                    resn = line[17:20].strip()
                    if resn == "PO4":
                        po4 += 1
                elif rec == "ENDMDL":
                    break
    except Exception:
        pass
    return mods, po4


def read_chain_seq(path, chain, max_lines=None):
    """按 ATOM 记录返回 (seq, [resseq], resnames)。按出现顺序，跳过非标准。"""
    seq, ids, names = [], [], []
    seen = set()
    try:
        with open(path, errors="ignore") as f:
            for line in f:
                if line.startswith("ENDMDL"):
                    break
                # 修饰残基（SEP/TPO/PTR 等）在 PDB 中以 HETATM 记录，必须一并纳入，
                # 否则磷酸化位点会从链序列中缺失，导致编号映射错位。
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                if line[21] != chain:
                    continue
                if line[17:20].strip() not in AA3 and line.startswith("HETATM"):
                    # 仅纳入氨基酸类修饰残基，排除配体/离子/水
                    if line[17:20].strip() not in ("SEP", "TPO", "PTR", "NEP", "MLY", "KCX"):
                        continue
                if line[12:16].strip() != "CA":
                    continue
                resn = line[17:20].strip()
                seqid = line[22:26].strip()
                ins = line[26]
                try:
                    n = int(seqid)
                except ValueError:
                    continue
                key = (n, ins)
                if key in seen:
                    continue
                seen.add(key)
                aa = AA3.get(resn)
                if aa is None:
                    # 磷酸化残基映射回母体氨基酸
                    aa = {"SEP": "S", "TPO": "T", "PTR": "Y", "NEP": "H"}.get(resn, "X")
                seq.append(aa); ids.append(n); names.append(resn)
    except Exception:
        pass
    return "".join(seq), ids, names


def identity(a, b):
    """全局简单一致性（等长时逐位比较；否则用滑动匹配兜底）"""
    if not a or not b:
        return 0.0, 0
    if len(a) == len(b):
        m = sum(1 for x, y in zip(a, b) if x == y)
        return m / len(a), m
    from difflib import SequenceMatcher
    sm = SequenceMatcher(None, a, b, autojunk=False)
    m = sum(b.size for b in sm.get_matching_blocks())
    return m / max(len(a), len(b)), m


# ---------- 主流程 ----------
def audit_protein(name, info, uni, args):
    gene = info["gene"]; acc = info["uniprot"]
    useq = uni.get(acc, {}).get("sequence", "")
    chain = info.get("chain", "A")
    pdb_chains = info.get("pdb_chains", {})
    sites = info.get("phos_sites", [])
    pdbs = info.get("pdbs", {})
    rec = {"protein": name, "gene": gene, "uniprot": acc, "chain": chain,
           "declared_sites": sites, "pdbs": OrderedDict(), "flags": [], "verdict": "PASS"}

    def flag(level, code, msg):
        rec["flags"].append({"level": level, "code": code, "msg": msg})
        if level == "FAIL":
            rec["verdict"] = "FAIL"
        elif level == "WARN" and rec["verdict"] == "PASS":
            rec["verdict"] = "WARN"

    for state, ids in pdbs.items():
        for pid in ids:
            f = PDBDIR / f"{pid.upper()}.pdb"
            e = {"pdb": pid.upper(), "state": state, "exists": f.exists()}
            if not f.exists():
                flag("FAIL", "C0_MISSING", f"{pid} 文件缺失")
                rec["pdbs"][pid] = e; continue

            title, mols = read_header(f)
            e["title"] = title[:110]
            e["molecules"] = mols

            # C1 身份核验：声明基因名/描述是否出现在头部
            blob = (title + " " + " ".join(mols)).upper()
            gene_ok = gene.upper() in blob
            desc_ok = any(w.upper() in blob for w in
                          re.findall(r"[A-Za-z]{4,}", info.get("description", ""))[:6])
            e["c1_identity"] = "HIT" if (gene_ok or desc_ok) else "MISS"
            if not (gene_ok or desc_ok):
                flag("FAIL", "C1_IDENTITY",
                     f"{pid} 头部未见 '{gene}' 或描述关键词；实际标示: "
                     f"{(mols[0] if mols else title[:60])}")

            # C2 真实磷酸化
            mods, po4 = read_modres(f)
            this_chain_c2 = pdb_chains.get(pid, pdb_chains.get(pid.upper(), chain))
            chain_mods = mods.get(this_chain_c2, {})
            other_mods = {c: v for c, v in mods.items() if c != chain}
            e["c2_modres_chain"] = {str(k): v[0] for k, v in chain_mods.items()}
            e["c2_modres_other_chains"] = {c: {str(k): v[0] for k, v in d.items()}
                                           for c, d in other_mods.items()}
            e["c2_po4_ligands"] = po4
            has_phos = bool(chain_mods) or bool(other_mods)
            e["c2_has_real_phos"] = has_phos
            if state == "phos" and not has_phos:
                flag("FAIL", "C2_NO_PHOS",
                     f"{pid} 声明为磷酸化态，但无 SEP/TPO/PTR 修饰记录"
                     + (f"（仅 {po4} 个 PO4 配体）" if po4 else ""))

            # C3 序列一致（每结构用各自链 ID）
            this_chain = pdb_chains.get(pid, pdb_chains.get(pid.upper(), chain))
            seq, ids_, names_ = read_chain_seq(f, this_chain)
            e["c3_chain"] = this_chain
            e["c3_chain_len"] = len(seq)
            e["c3_range"] = f"{min(ids_)}..{max(ids_)}" if ids_ else "-"
            if not seq:
                flag("FAIL", "C3_NO_CHAIN", f"{pid} 链 {this_chain} 无 ATOM 记录")
                rec["pdbs"][pid] = e; continue
            cov, nmatch = identity(seq, useq)
            e["c3_cov_vs_uniprot"] = round(cov, 3)
            e["c3_nmatch"] = nmatch
            if seq not in useq and cov < args.min_cov:
                flag("FAIL", "C3_SEQ_MISMATCH",
                     f"{pid} 链 {this_chain} 与 UniProt 覆盖仅 {cov:.0%}（阈值 {args.min_cov:.0%}）")

            # C4 覆盖范围（供配对重叠判定）
            e["_ids"] = (min(ids_), max(ids_)) if ids_ else None
            e["_seq"] = seq
            rec["pdbs"][pid] = e

    # ---- C4：配对重叠判定（仅对有两态 PDB 的蛋白） ----
    ent = list(rec["pdbs"].items())
    ph = [v for k, v in ent if v["state"] == "phos" and v.get("_ids")]
    up = [v for k, v in ent if v["state"] == "unphos" and v.get("_ids")]
    rec["c4"] = {}
    if ph and up:
        a, b = ph[0], up[0]
        lo = max(a["_ids"][0], b["_ids"][0]); hi = min(a["_ids"][1], b["_ids"][1])
        ov = max(0, hi - lo + 1)
        rec["c4"] = {"phos_range": list(a["_ids"]), "unphos_range": list(b["_ids"]),
                     "overlap_residues": ov}
        if ov < args.min_overlap:
            flag("FAIL", "C4_NO_OVERLAP",
                 f"phos({a['pdb']} {a['_ids'][0]}-{a['_ids'][1]}) 与 "
                 f"unphos({b['pdb']} {b['_ids'][0]}-{b['_ids'][1]}) 重叠仅 {ov} 残基")
        # 位点在重叠区内？
        site_nums = [int("".join(c for c in s if c.isdigit())) for s in sites if any(c.isdigit() for c in s)]
        in_ph = [n for n in site_nums if a["_ids"][0] <= n <= a["_ids"][1]]
        in_up = [n for n in site_nums if b["_ids"][0] <= n <= b["_ids"][1]]
        rec["c4"]["sites_in_phos_range"] = in_ph
        rec["c4"]["sites_in_unphos_range"] = in_up
        if site_nums and not in_ph:
            flag("WARN", "C4_SITE_NOT_IN_PHOS",
                 f"声明位点 {sites} 不在 phos 结构 {a['pdb']} 覆盖范围内（编号体系可能不同）")

    # ---- C5：配对纯度（两态结构的链组成是否一致）----
    if ph and up:
        try:
            f_up = PDBDIR / f"{up[0]['pdb']}.pdb"
            f_ph = PDBDIR / f"{ph[0]['pdb']}.pdb"
            ch_up = {ln[21] for ln in open(f_up, errors="ignore") if ln.startswith("ATOM")}
            ch_ph = {ln[21] for ln in open(f_ph, errors="ignore") if ln.startswith("ATOM")}
            rec["c5_chains_unphos"] = sorted(ch_up)
            rec["c5_chains_phos"] = sorted(ch_ph)
            if ch_up != ch_ph:
                flag("WARN", "C5_COMPOSITION",
                     f"两态链组成不同 unphos={sorted(ch_up)} vs phos={sorted(ch_ph)}"
                     f"（复合物/构建体差异可能混入磷酸化效应）")
            rec["c5_pure_pair"] = (ch_up == ch_ph)
        except Exception:
            rec["c5_pure_pair"] = None

    # 清理内部字段
    for v in rec["pdbs"].values():
        v.pop("_ids", None); v.pop("_seq", None)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cov", type=float, default=0.5, help="与 UniProt 覆盖度阈值")
    ap.add_argument("--min-overlap", type=int, default=20, help="两态最少重叠残基数")
    args = ap.parse_args()

    uni_p = _ROOT / "data" / "curated" / "uniprot_sequences.json"
    uni = json.loads(uni_p.read_text()) if uni_p.exists() else {}

    results = [audit_protein(n, i, uni, args) for n, i in BENCHMARK_PROTEINS.items()]

    # ---- 写 JSON ----
    (OUT/"pair_audit.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))

    # ---- 控制台汇总 ----
    print("="*104)
    print("  基准配对质量筛查")
    print("="*104)
    print(f"{'蛋白':<24}{'判定':<7}{'phos PDB':<10}{'真实磷酸化':<11}{'unphos PDB':<11}{'重叠':<7}主要问题")
    print("-"*104)
    for r in results:
        ph = next((v for v in r["pdbs"].values() if v["state"]=="phos"), {})
        up = next((v for v in r["pdbs"].values() if v["state"]=="unphos"), {})
        phos_ok = "是" if ph.get("c2_has_real_phos") else ("-" if not ph else "否")
        ov = r["c4"].get("overlap_residues", "-") if r.get("c4") else "-"
        fails = [f["code"] for f in r["flags"] if f["level"]=="FAIL"]
        print(f"{r['protein']:<24}{r['verdict']:<7}{ph.get('pdb','-'):<10}"
              f"{phos_ok:<11}{up.get('pdb','-'):<11}{str(ov):<7}{','.join(fails[:3])}")

    n_pass = sum(1 for r in results if r["verdict"]=="PASS")
    n_warn = sum(1 for r in results if r["verdict"]=="WARN")
    n_fail = sum(1 for r in results if r["verdict"]=="FAIL")
    print("-"*104)
    print(f"  合计 {len(results)} 蛋白：PASS {n_pass} | WARN {n_warn} | FAIL {n_fail}")

    # ---- 写 Markdown ----
    md = ["# 基准配对质量筛查报告", "",
          f"筛选阈值：与 UniProt 覆盖度 ≥ {args.min_cov:.0%}；两态重叠 ≥ {args.min_overlap} 残基。", "",
          "## 汇总", "",
          "| 蛋白 | 判定 | phos PDB | 真实磷酸化 | unphos PDB | 重叠残基 | 问题 |",
          "|---|---|---|---|---|---|---|"]
    for r in results:
        ph = next((v for v in r["pdbs"].values() if v["state"]=="phos"), {})
        up = next((v for v in r["pdbs"].values() if v["state"]=="unphos"), {})
        phos_ok = "是" if ph.get("c2_has_real_phos") else ("—" if not ph else "**否**")
        ov = r["c4"].get("overlap_residues", "—") if r.get("c4") else "—"
        probs = "; ".join(f"{f['code']}" for f in r["flags"]) or "—"
        md.append(f"| {r['protein']} | {r['verdict']} | {ph.get('pdb','—')} | {phos_ok} "
                  f"| {up.get('pdb','—')} | {ov} | {probs} |")
    md += ["", "## 逐条问题明细", ""]
    for r in results:
        if not r["flags"]: continue
        md.append(f"### {r['protein']}（{r['verdict']}）")
        for fl in r["flags"]:
            md.append(f"- **[{fl['level']}] {fl['code']}** — {fl['msg']}")
        md.append("")
    (OUT/"pair_audit.md").write_text("\n".join(md))
    print(f"\n[OK] {OUT}/pair_audit.json")
    print(f"[OK] {OUT}/pair_audit.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
