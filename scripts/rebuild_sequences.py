#!/usr/bin/env python3
"""
rebuild_sequences.py — P0-1 序列重建：以 UniProt 为唯一权威序列源
====================================================================
问题背景（见审查报告 2.1 节）：
  原始 data/fasta/*.fasta 由 BioPython 从 PDB **观测坐标**提取（seq1() + 标准残基过滤），
  导致 (a) 长度系统性偏短、(b) phos/unphos 两态不可比（IRF3 239 vs 110 aa）、
  (c) 无法与 UniProt 编号的磷酸化位点对齐。

本脚本：
  1. 从 UniProt REST 拉取 18 个基准蛋白的全长 canonical 序列
  2. 校验 benchmark_dataset.py 中每个 phos_site 在 UniProt 序列上的残基类型（S/T/Y）
  3. 解析本地 PDB 的 DBREF/SEQADV，建立 PDB 编号 <-> UniProt 编号映射
  4. 写出 data/curated/{gene}_{state}.fasta（UniProt 全长）与
     data/curated/uniprot_sequences.json、data/curated/pdb_uniprot_mapping.json

用法：
  python scripts/rebuild_sequences.py            # 全部重建
  python scripts/rebuild_sequences.py --verify   # 仅校验，不写文件
"""
import sys, json, argparse, urllib.request, time
from pathlib import Path
from collections import OrderedDict

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

CURATED = _ROOT / "data" / "curated"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{acc}.json"
PTM_AA = {"SEP": "S", "TPO": "T", "PTR": "Y"}


def fetch_uniprot_seq(acc, retries=3):
    """拉取 UniProt canonical FASTA 序列。"""
    for i in range(retries):
        try:
            with urllib.request.urlopen(UNIPROT_FASTA.format(acc=acc), timeout=20) as r:
                txt = r.read().decode()
            lines = txt.splitlines()
            return "".join(l.strip() for l in lines[1:]), lines[0]
        except Exception as e:
            if i == retries - 1:
                print(f"    [FAIL] {acc}: {type(e).__name__} {e}")
                return None, None
            time.sleep(2)
    return None, None


def map_pdb_to_uniprot(pdb_id, acc, chain):
    """
    用 PDB 的 DBREF 记录建立 PDB 残基编号 -> UniProt 编号映射。
    返回 (offset, evidence) ；offset = uniprot_pos - pdb_pos。
    """
    from Bio.PDB import PDBParser
    pdb_file = _ROOT / "data" / "pdb" / f"{pdb_id.upper()}.pdb"
    if not pdb_file.exists():
        return None, "PDB file missing"
    offset = None
    evidence = []
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("DBREF") and line[12:16].strip() == chain:
                    db_acc = line[33:41].strip()
                    if db_acc.upper() == acc.upper():
                        pdb_start = int(line[20:25])
                        seq_start = int(line[42:49])
                        offset = seq_start - pdb_start
                        evidence.append(f"DBREF {pdb_id}:{chain} pdb{pdb_start}->uniprot{seq_start}")
                elif line.startswith("SEQADV") and line[12:16].strip() == chain:
                    db_acc = line[21:29].strip()
                    if db_acc.upper() == acc.upper():
                        evidence.append(f"SEQADV {line[18:22].strip()} {line[42:49].strip()}")
    except Exception as e:
        return None, f"parse error: {e}"
    return offset, "; ".join(evidence) if evidence else "no DBREF/SEQADV for this chain"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true", help="仅校验不写文件")
    args = ap.parse_args()
    CURATED.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("  P0-1 序列重建：UniProt canonical 为唯一权威源")
    print("=" * 72)

    seqs = OrderedDict()
    reports = OrderedDict()
    n_ok = n_bad = 0

    for name, info in BENCHMARK_PROTEINS.items():
        acc = info["uniprot"]
        gene = info["gene"]
        print(f"\n[{name}] {gene} / {acc}")
        seq, header = fetch_uniprot_seq(acc)
        if not seq:
            n_bad += 1
            reports[name] = {"status": "FETCH_FAILED", "uniprot": acc}
            continue
        n_ok += 1
        seqs[acc] = {"protein": name, "gene": gene, "sequence": seq,
                     "length": len(seq), "uniprot_header": header}

        # --- 校验磷酸化位点 ---
        site_checks = []
        for site in info.get("phos_sites", []):
            aa = site[0]
            try:
                pos = int("".join(c for c in site if c.isdigit()))
            except ValueError:
                site_checks.append({"site": site, "status": "UNPARSEABLE"})
                continue
            if pos < 1 or pos > len(seq):
                site_checks.append({"site": site, "pos": pos,
                                    "status": f"OUT_OF_RANGE (len={len(seq)})"})
                continue
            actual = seq[pos - 1]
            ok = (actual == aa)
            site_checks.append({"site": site, "pos": pos, "expected": aa,
                                "actual": actual, "status": "OK" if ok else "MISMATCH"})
            flag = "OK " if ok else "!!!"
            print(f"   {flag} {site:<10s} UniProt[{pos}] = {actual} (expect {aa})")

        # --- PDB 编号映射 ---
        pdbs = info.get("pdbs", {})
        chain = info.get("chain", "A")
        mappings = {}
        for state, ids in pdbs.items():
            for pid in ids:
                off, ev = map_pdb_to_uniprot(pid, acc, chain)
                mappings[pid] = {"state": state, "chain": chain,
                                 "offset_uniprot_minus_pdb": off, "evidence": ev}
                if off is not None:
                    print(f"   [map] {pid} ({state}): offset {off:+d} | {ev[:60]}")

        reports[name] = {
            "status": "OK", "gene": gene, "uniprot": acc,
            "uniprot_length": len(seq),
            "declared_length_in_benchmark": info.get("length"),
            "length_match": (info.get("length") == len(seq)),
            "sites": site_checks,
            "pdb_mappings": mappings,
        }

    print("\n" + "=" * 72)
    print(f"  拉取结果：成功 {n_ok} / 失败 {n_bad}")
    print("=" * 72)

    if args.verify:
        print("\n[--verify] 未写文件。")
        return

    # 写 UniProt 序列库
    (CURATED / "uniprot_sequences.json").write_text(
        json.dumps(seqs, indent=2, ensure_ascii=False))

    # 写逐蛋白报告
    (CURATED / "sequence_audit.json").write_text(
        json.dumps(reports, indent=2, ensure_ascii=False))

    # 写 FASTA（每蛋白 unphos + phos，序列同为 UniProt 全长；差异体现在修饰注释）
    n_fasta = 0
    for name, rep in reports.items():
        if rep["status"] != "OK":
            continue
        gene = rep["gene"]
        seq = seqs[rep["uniprot"]]["sequence"]
        for state in ("unphos", "phos"):
            out = CURATED / f"{gene}_{state}.fasta"
            hdr = f">{gene}_{state}|{rep['uniprot']}|uniprot_len={len(seq)}"
            if state == "phos":
                sites = ",".join(f"{s['site']}" for s in rep["sites"] if s["status"] == "OK")
                hdr += f"|ptm_sites={sites}"
            out.write_text(f"{hdr}\n{seq}\n")
            n_fasta += 1
    print(f"\n[OK] 写出 {n_fasta} 条 FASTA 到 {CURATED}/")

    # 汇总映射表
    all_map = {p: r.get("pdb_mappings", {}) for p, r in reports.items() if r["status"] == "OK"}
    (CURATED / "pdb_uniprot_mapping.json").write_text(
        json.dumps(all_map, indent=2, ensure_ascii=False))
    print(f"[OK] 写出 PDB->UniProt 映射到 {CURATED}/pdb_uniprot_mapping.json")

    # 打印长度核对表
    print("\n" + "-" * 72)
    print(f"{'protein':<22s}{'benchmark_len':>14s}{'uniprot_len':>13s}  match")
    print("-" * 72)
    for name, rep in reports.items():
        if rep["status"] != "OK":
            print(f"{name:<22s}{'--':>14s}{'FETCH FAILED':>13s}")
            continue
        decl = rep.get("declared_length_in_benchmark")
        m = "OK" if rep["length_match"] else (f"DIFF {decl}" if decl else "no-declared")
        print(f"{name:<22s}{str(decl if decl is not None else '-'):>14}"
              f"{rep['uniprot_length']:>13}  {m}")


if __name__ == "__main__":
    main()
