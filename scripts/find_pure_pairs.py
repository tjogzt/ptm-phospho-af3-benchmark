#!/usr/bin/env python3
"""
find_pure_pairs.py — 搜索"复合物组成完全一致"的磷酸化配对（C5 纯净配对）
==========================================================================
动机：10 组配对中仅 2 组纯净，其余混入复合物/配体差异，严重稀释磷酸化信号。

两阶段策略（避免逐个下载 PDB）：
  Stage 1  Search API 列出某 UniProt 的全部 entry
           → GraphQL 批量取每个 entry 的 **UniProt ID 集合** 与分辨率
           → 以 UniProt 集合为"复合物组成"键分组（比分子名可靠）
  Stage 2  仅在"组内既有 phos 又有 unphos 候选"的组里，下载 PDB 头部
           解析 MODRES 判定磷酸化状态与位点

输出：audit/pure_pair_candidates.json

用法：
  python3 scripts/find_pure_pairs.py --all
  python3 scripts/find_pure_pairs.py --only CDK2 CDK1 AMPK_alpha2
"""
import sys, json, argparse, re, time
from pathlib import Path
from collections import defaultdict
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

OUT = _ROOT / "audit"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
FILES = "https://files.rcsb.org/download/{pid}.pdb"
GRAPHQL = "https://data.rcsb.org/graphql"
PHOS = {"SEP", "TPO", "PTR", "NEP"}
HEADER_BYTES = 80000


def search_entries(acc, rows=200, retries=3):
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers"
                         ".reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows},
                             "results_content_type": ["experimental"]}}
    for attempt in range(retries):
        req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(),
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=40) as r:
                got = [x["identifier"] for x in
                       json.loads(r.read().decode()).get("result_set", [])]
            if got:
                return got
        except Exception:
            pass
        time.sleep(3 * (attempt + 1))          # 指数退避
    return []


def gql_composition(pids):
    """批量取 UniProt 组成 + 分辨率"""
    if not pids: return {}
    ids = json.dumps(pids).replace('"', '"')
    q = ('{ entries(entry_ids: %s) { rcsb_id rcsb_entry_info { resolution_combined } '
         'polymer_entities { rcsb_polymer_entity_container_identifiers { '
         'reference_sequence_identifiers { database_accession database_name } } } } }') % ids
    req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": q}).encode(),
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=45) as r:
            d = json.loads(r.read().decode())
    except Exception:
        return {}
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
        out[e["rcsb_id"]] = {"uniprot_set": frozenset(unis),
                             "resolution": res[0] if res else None}
    return out


def fetch_header(pid):
    req = urllib.request.Request(FILES.format(pid=pid),
                                 headers={"Range": f"bytes=0-{HEADER_BYTES}"})
    try:
        with urllib.request.urlopen(req, timeout=25) as r:
            return r.read().decode("utf-8", "ignore")
    except Exception:
        return ""


def parse_modres(txt):
    title, mods = [], []
    for line in txt.splitlines():
        if line.startswith("TITLE"): title.append(line[10:].strip())
        elif line.startswith("MODRES"):
            if line[12:15].strip() in PHOS:
                mods.append(int(line[18:22].strip() or 0))
    return " ".join(title)[:95], sorted(set(mods))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--max", type=int, default=120)
    args = ap.parse_args()

    targets = ([k for k, v in BENCHMARK_PROTEINS.items() if v.get("paired")]
               if args.all or not args.only else args.only)
    report = {}
    print("="*110)
    print("  纯净配对搜索（组成键 = UniProt ID 集合）")
    print("="*110)

    for name in targets:
        info = BENCHMARK_PROTEINS.get(name)
        if not info: continue
        acc = info["uniprot"]
        sites = [int("".join(c for c in s if c.isdigit())) for s in info.get("phos_sites", [])]
        entries = search_entries(acc)[:args.max]
        comp = {}
        for i in range(0, len(entries), 60):
            comp.update(gql_composition(entries[i:i+60]))
        # 组内需有多于 1 个成员才可能配对
        groups = defaultdict(list)
        for pid, c in comp.items():
            if c["uniprot_set"]:
                groups[c["uniprot_set"]].append((pid, c["resolution"]))
        cand_groups = {k: v for k, v in groups.items() if len(v) > 1}

        print(f"\n### {name} ({acc}) 目标位点 {sites} | entry {len(entries)} | "
              f"组成组 {len(groups)}（可配对 {len(cand_groups)}）")

        # Stage 2：只查候选组
        to_check = [pid for v in cand_groups.values() for pid, _ in v]
        phos_state = {}
        if to_check:
            with ThreadPoolExecutor(max_workers=10) as ex:
                futs = {ex.submit(fetch_header, pid): pid for pid in to_check}
                for fu in as_completed(futs):
                    pid = futs[fu]
                    txt = fu.result()
                    if not txt: continue
                    title, mods = parse_modres(txt)
                    phos_state[pid] = {"sites": mods, "title": title,
                                       "res": comp.get(pid, {}).get("resolution")}

        pairs = []
        for ukey, members in groups.items():
            ph = [(p, r) for p, r in members if phos_state.get(p, {}).get("sites")]
            un = [(p, r) for p, r in members
                  if p in phos_state and not phos_state[p]["sites"]]
            if not (ph and un): continue
            hit = any(abs(s-t) <= 2 for p, _ in ph
                      for s in phos_state[p]["sites"] for t in sites)
            pairs.append({
                "composition": sorted(ukey), "n_components": len(ukey),
                "phos": sorted([{"pdb": p, "sites": phos_state[p]["sites"],
                                 "res": r} for p, r in ph], key=lambda x: x["res"] or 99),
                "unphos": sorted([{"pdb": p, "res": r} for p, r in un],
                                 key=lambda x: x["res"] or 99),
                "target_site_hit": hit,
            })
        pairs.sort(key=lambda x: (not x["target_site_hit"], x["n_components"]))
        report[name] = {"uniprot": acc, "sites": sites,
                        "n_entries": len(entries), "pairs": pairs}
        if not pairs:
            print("    ✗ 无同组成配对")
        for x in pairs[:3]:
            tag = " ★位点匹配" if x["target_site_hit"] else ""
            print(f"    ✓ {x['n_components']} 组分 {x['composition'][:3]}{tag}")
            print(f"        phos  : {[(p['pdb'], p['res'], p['sites']) for p in x['phos'][:3]]}")
            print(f"        unphos: {[(u['pdb'], u['res']) for u in x['unphos'][:3]]}")

    outp = OUT / "pure_pair_candidates.json"
    merged = {}
    if outp.exists():
        try: merged = json.loads(outp.read_text())
        except Exception: pass
    for k, v in report.items():
        # frozenset 不能直接存 JSON，composition 已转 list
        merged[k] = v
    outp.write_text(json.dumps(merged, indent=2, ensure_ascii=False))
    print(f"\n[OK] {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
