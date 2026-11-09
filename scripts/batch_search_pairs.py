#!/usr/bin/env python3
"""
batch_search_pairs.py — 批量筛选"链组成一致"的磷酸化纯净配对
===============================================================
输入：候选蛋白 (UniProt accession, 激活环磷酸化位点) 列表
输出：每个候选是否存在纯净配对，及推荐结构

流程（同 find_pure_pairs.py，但支持任意 UniProt）：
  1. Search API 列出 entry
  2. GraphQL 批量取 UniProt 组成 + 分辨率
  3. 按组成分组，组内含 phos 与 unphos → 候选
  4. 下载头部确认磷酸化状态与位点

用法：python3 scripts/batch_search_pairs.py --batch 1
"""
import sys, json, argparse, re, time
from pathlib import Path
from collections import defaultdict
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

_ROOT = Path(__file__).resolve().parent.parent
OUT = _ROOT / "audit"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
FILES = "https://files.rcsb.org/download/{pid}.pdb"
GRAPHQL = "https://data.rcsb.org/graphql"
PHOS = {"SEP", "TPO", "PTR", "NEP"}
HEADER_BYTES = 80000

# 候选：经典磷酸化调控蛋白（激酶激活环为主）
CANDIDATES = [
    ("MAPK3",   "P27361", [202, 204]),   # ERK1
    ("MAPK8",   "P45983", [183, 185]),   # JNK1
    ("MAPK9",   "P45984", [183, 185]),   # JNK2
    ("MAPK11",  "Q15759", [180, 182]),   # p38β
    ("GSK3B",   "P49841", [216]),        # GSK3β Y216
    ("CHEK1",   "O14757", [317, 345]),   # Chk1
    ("AURKB",   "Q96GD4", [232]),        # Aurora B
    ("PLK1",    "P53350", [210]),        # PLK1
    ("PAK1",    "Q13153", [423]),        # PAK1
    ("ABL1",    "P00519", [245]),        # ABL1
    ("CDK5",    "Q00535", [15]),         # CDK5 (p35 激活)
    ("PRKAA1",  "Q13131", [183]),        # AMPK α1
    ("RPS6KA1", "Q15418", [359, 363, 380]),  # p90RSK
    ("CHEK2",   "O96017", [68]),         # Chk2
    ("WEE1",    "P30291", [53]),         # Wee1
    ("AKT1",    "P31749", [308, 473]),   # AKT1
    ("PDPK1",   "O15530", [9]),          # PDK1
    ("TBK1",    "Q9UHD2", [172]),        # TBK1
    ("STK11",   "Q15831", [189]),        # LKB1
    ("IKBKB",   "O14920", [177, 181]),   # IKKβ
]


def search_entries(acc, rows=120, retries=2):
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers"
                         ".reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows},
                             "results_content_type": ["experimental"]}}
    for a in range(retries):
        req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(),
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=35) as r:
                got = [x["identifier"] for x in json.loads(r.read().decode()).get("result_set", [])]
            if got: return got
        except Exception: pass
        time.sleep(2*(a+1))
    return []


def gql_comp(pids):
    if not pids: return {}
    q = ('{ entries(entry_ids: %s) { rcsb_id rcsb_entry_info { resolution_combined } '
         'polymer_entities { rcsb_polymer_entity_container_identifiers { '
         'reference_sequence_identifiers { database_accession database_name } } } } }') % json.dumps(pids)
    req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": q}).encode(),
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=45) as r:
            d = json.loads(r.read().decode())
    except Exception: return {}
    out = {}
    for e in d.get("data", {}).get("entries") or []:
        if not e: continue
        unis = set()
        for pe in e.get("polymer_entities") or []:
            for rid in (pe["rcsb_polymer_entity_container_identifiers"]
                        .get("reference_sequence_identifiers") or []):
                if rid.get("database_name") == "UniProt":
                    unis.add(rid["database_accession"])
        res = (e.get("rcsb_entry_info") or {}).get("resolution_combined")
        out[e["rcsb_id"]] = {"set": frozenset(unis), "res": res[0] if res else None}
    return out


def fetch_modres(pid):
    req = urllib.request.Request(FILES.format(pid=pid),
                                 headers={"Range": f"bytes=0-{HEADER_BYTES}"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            txt = r.read().decode("utf-8", "ignore")
    except Exception: return None
    mods = []
    for line in txt.splitlines():
        if line.startswith("MODRES") and line[12:15].strip() in PHOS:
            mods.append(int(line[18:22].strip() or 0))
    return sorted(set(mods))


def probe(name, acc, sites):
    entries = search_entries(acc, rows=100)
    if len(entries) < 2: return None
    comp = {}
    for i in range(0, len(entries), 60):
        comp.update(gql_comp(entries[i:i+60]))
    groups = defaultdict(list)
    for pid, c in comp.items():
        if c["set"]: groups[c["set"]].append((pid, c["res"]))
    cand = {k: v for k, v in groups.items() if len(v) > 1}
    if not cand: return None
    # 只检查候选组内的结构
    to_check = [p for v in cand.values() for p, _ in v]
    state = {}
    with ThreadPoolExecutor(max_workers=10) as ex:
        futs = {ex.submit(fetch_modres, p): p for p in to_check}
        for fu in as_completed(futs):
            p = futs[fu]; m = fu.result()
            if m is not None: state[p] = m
    best = None
    for ukey, members in groups.items():
        ph = [(p, r) for p, r in members if state.get(p)]
        un = [(p, r) for p, r in members if p in state and not state[p]]
        if not (ph and un): continue
        # 位点匹配筛选
        hit = [(p, r, state[p]) for p, r in ph
               if any(abs(s-t) <= 2 for s in state[p] for t in sites)]
        if not hit: continue
        hit.sort(key=lambda x: x[1] or 99)
        un.sort(key=lambda x: x[1] or 99)
        rec = {"composition": sorted(ukey), "n": len(ukey),
               "phos": [{"pdb": p, "res": r, "sites": m} for p, r, m in hit[:3]],
               "unphos": [{"pdb": p, "res": r} for p, r in un[:3]]}
        if best is None or len(ukey) < best["n"]: best = rec
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=int, default=0, help="0=全部, 1..N=第N批(每批5个)")
    args = ap.parse_args()
    cands = CANDIDATES
    if args.batch > 0:
        cands = CANDIDATES[(args.batch-1)*5: args.batch*5]

    report = {}
    print("="*104)
    print(f"  批量纯净配对搜索（{len(cands)} 个候选）")
    print("="*104)
    for name, acc, sites in cands:
        try:
            r = probe(name, acc, sites)
        except Exception as e:
            print(f"  {name:<10} ERROR {type(e).__name__}"); continue
        if r:
            print(f"  ✓ {name:<10} ({acc}) 位点{sites} → 组分 {r['n']} | "
                  f"phos={[p['pdb'] for p in r['phos']]} unphos={[u['pdb'] for u in r['unphos']]}")
            report[name] = {"acc": acc, "sites": sites, **r}
        else:
            print(f"  ✗ {name:<10} ({acc}) 无纯净配对")
        time.sleep(0.5)

    outp = OUT/"batch_pure_candidates.json"
    merged = {}
    if outp.exists():
        try: merged = json.loads(outp.read_text())
        except Exception: pass
    merged.update(report)
    outp.write_text(json.dumps(merged, indent=2, ensure_ascii=False))
    print(f"\n[OK] {outp}（累计 {len(merged)} 个有纯净配对的候选）")
    return report


if __name__ == "__main__":
    main()
