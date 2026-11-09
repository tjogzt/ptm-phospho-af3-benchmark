#!/usr/bin/env python3
"""
find_new_pairs.py — 为扩充基准自动筛找"磷酸化/非磷酸化"配对完整的候选蛋白
============================================================================
对候选 UniProt 逐一：
  1. 检索全部 PDB entry
  2. 按是否含 SEP/TPO/PTR 分成 phos / unphos 两组
  3. 检查是否覆盖该蛋白的经典激活环位点
  4. 输出可配对者（两侧均存在且有覆盖）

用法：python3 scripts/find_new_pairs.py --acc P06493 --sites 161
"""
import sys, json, argparse, re
from pathlib import Path
import urllib.request
from difflib import SequenceMatcher

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "scripts"))
from audit_pairs import read_header, read_modres, read_chain_seq

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
FILES = "https://files.rcsb.org/download/{pid}.pdb"


def search(acc, rows=120):
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers"
                         ".reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows},
                             "results_content_type": ["experimental"]}}
    req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=40) as r:
        d = json.loads(r.read().decode())
    return [x["identifier"] for x in d.get("result_set", [])], d.get("total_count", 0)


def get_local(pid):
    f = _ROOT / "data" / "pdb" / f"{pid}.pdb"
    if f.exists() and f.stat().st_size > 20000: return f
    try:
        with urllib.request.urlopen(FILES.format(pid=pid), timeout=45) as r:
            data = r.read()
        if len(data) < 10000: return None
        f.write_bytes(data); return f
    except Exception:
        return None


def uniprot_seq(acc):
    with urllib.request.urlopen(
        f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=25) as r:
        lines = r.read().decode().splitlines()
    return "".join(l.strip() for l in lines[1:])


def probe(pid, useq, sites):
    f = get_local(pid)
    if not f: return None
    mods, _ = read_modres(f)
    phos = sorted({s for c, d in mods.items() for s in d})
    best = None
    for ch in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        try: seq, ids, _ = read_chain_seq(f, ch)
        except Exception: continue
        if len(seq) < 50: continue
        sm = SequenceMatcher(None, seq, useq, autojunk=False)
        nm = sum(b.size for b in sm.get_matching_blocks())
        cov = nm / max(len(seq), len(useq))
        if best is None or nm > best[3]:
            best = (ch, (min(ids), max(ids)) if ids else None, cov, nm)
    if not best: return None
    ch, rng, cov, _ = best
    if cov < 0.35: return None
    covered = [s for s in sites if rng and rng[0] <= s <= rng[1]]
    if not covered: return None
    title, _ = read_header(f)
    return {"pdb": pid, "chain": ch, "range": list(rng), "cov": round(cov, 3),
            "phos": phos, "sites_covered": covered, "title": title[:95]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--acc", required=True)
    ap.add_argument("--sites", nargs="+", type=int, required=True)
    ap.add_argument("--max", type=int, default=60)
    args = ap.parse_args()

    useq = uniprot_seq(args.acc)
    entries, total = search(args.acc)
    print(f"UniProt {args.acc} len={len(useq)} | entry {total} | 检查 {min(len(entries), args.max)}")
    ph, un = [], []
    for pid in entries[:args.max]:
        r = probe(pid, useq, args.sites)
        if not r: continue
        (ph if r["phos"] else un).append(r)
    print(f"\n  含磷酸化: {len(ph)} 个")
    for r in sorted(ph, key=lambda x: x["pdb"])[:5]:
        print(f"    {r['pdb']} 链{r['chain']} cov={r['cov']:.0%} sites={r['phos']} | {r['title'][:60]}")
    print(f"  无磷酸化: {len(un)} 个")
    for r in sorted(un, key=lambda x: x["pdb"])[:5]:
        print(f"    {r['pdb']} 链{r['chain']} cov={r['cov']:.0%} | {r['title'][:60]}")
    print(f"\n  可配对: {'是' if ph and un else '否'}")
    return {"acc": args.acc, "sites": args.sites, "phos": ph, "unphos": un}


if __name__ == "__main__":
    main()
