#!/usr/bin/env python3
"""
find_unphos_structures.py — 为非磷酸化态检索实验结构（补齐配对）
==================================================================
用途：为当前只有 phos 实验结构的蛋白，寻找**不含任何磷酸化修饰**的对照结构，
使其升级为成对结构，从而参与 AF3-vs-实验 位移验证。

筛选条件（四条）：
  U1 无 SEP/TPO/PTR 修饰（MODRES + 残基名双重检查）
  U2 结构覆盖目标磷酸化位点（否则无法比较该位点位移）
  U3 与 UniProt 序列覆盖度达标
  U4 分辨率记录（供人工裁决，不硬性过滤）

输出：audit/unphos_candidates.json + 控制台候选表

用法：
  python3 scripts/find_unphos_structures.py --only Glycogen_Phosphorylase
  python3 scripts/find_unphos_structures.py --only c-Src_kinase AMPK_alpha2
"""
import sys, json, argparse, time, re
from pathlib import Path
import urllib.request

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "data"))
sys.path.insert(0, str(_ROOT / "scripts"))
from benchmark_dataset import BENCHMARK_PROTEINS
from audit_pairs import read_header, read_modres, read_chain_seq

OUT = _ROOT / "audit"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
FILES = "https://files.rcsb.org/download/{pid}.pdb"
HEADER_BYTES = 262144
PHOS = {"SEP", "TPO", "PTR", "NEP"}


def search_entries(acc, rows=300):
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


def ensure_local(pid):
    """确保 PDB 文件在本地（只取头部用于 MODRES，但结构比较需完整文件）"""
    f = _ROOT / "data" / "pdb" / f"{pid}.pdb"
    if f.exists() and f.stat().st_size > 20000:
        return f, True
    try:
        req = urllib.request.Request(FILES.format(pid=pid))
        with urllib.request.urlopen(req, timeout=45) as r:
            data = r.read()
        if len(data) < 10000:
            return None, False
        f.write_bytes(data)
        return f, True
    except Exception:
        return None, False


def resolution_of(path):
    try:
        with open(path, errors="ignore") as fh:
            for line in fh:
                if line.startswith("REMARK   2 RESOLUTION"):
                    m = re.search(r"(\d+\.\d+)\s+ANGSTROM", line)
                    if m: return float(m.group(1))
                if line.startswith("ATOM"): break
    except Exception:
        pass
    return None


def evaluate(pid, acc, sites, useq):
    """返回评估结果 dict；不可用返回 None"""
    f, ok = ensure_local(pid)
    if not ok: return None
    mods, _po4 = read_modres(f)
    phos = {s for c, d in mods.items() for s in d}
    # U1: 无磷酸化
    if phos:
        return {"pdb": pid, "usable": False, "reason": f"含磷酸化 {sorted(phos)}"}
    # 找与 UniProt 最匹配的链
    best = None
    for ch in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789":
        try:
            seq, ids, _ = read_chain_seq(f, ch)
        except Exception:
            continue
        if len(seq) < 50: continue
        from difflib import SequenceMatcher
        sm = SequenceMatcher(None, seq, useq, autojunk=False)
        nm = sum(b.size for b in sm.get_matching_blocks())
        cov = nm / max(len(seq), len(useq))
        if best is None or nm > best[4]:
            best = (ch, len(seq), (min(ids), max(ids)) if ids else None, nm, cov)
    if not best:
        return {"pdb": pid, "usable": False, "reason": "无可用链"}
    ch, n, rng, nm, cov = best
    if cov < 0.3:
        return {"pdb": pid, "usable": False, "reason": f"序列覆盖仅 {cov:.0%}"}
    # U2: 覆盖目标位点
    sites_in = [s for s in sites if rng and rng[0] <= s <= rng[1]]
    title, _ = read_header(f)
    return {"pdb": pid, "usable": bool(sites_in), "chain": ch, "n_res": n,
            "range": list(rng) if rng else None, "cov": round(cov, 3),
            "sites_covered": sites_in, "sites_missing": [s for s in sites if s not in sites_in],
            "resolution": resolution_of(f), "title": title[:110],
            "reason": "" if sites_in else "未覆盖目标位点"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--max-entries", type=int, default=40, help="每蛋白最多检查的 entry 数")
    args = ap.parse_args()

    TARGETS = ["Glycogen_Phosphorylase", "c-Src_kinase", "AMPK_alpha2",
               "PKA_Calpha", "CTD_kinase"]
    targets = args.only or TARGETS
    uni = json.loads((_ROOT/"data"/"curated"/"uniprot_sequences.json").read_text())

    report = {}
    for name in targets:
        info = BENCHMARK_PROTEINS.get(name)
        if not info: 
            print(f"[SKIP] {name}"); continue
        acc = info["uniprot"]; useq = uni[acc]["sequence"]
        sites = [int("".join(c for c in s if c.isdigit())) for s in info.get("phos_sites", [])]
        print(f"\n### {name} ({acc}) 目标位点 {sites}")
        try:
            entries, total = search_entries(acc)
        except Exception as e:
            print(f"    检索失败: {e}"); continue
        print(f"    关联 entry {total} 个，检查前 {min(len(entries), args.max_entries)} 个")

        cands = []
        for pid in entries[:args.max_entries]:
            try:
                r = evaluate(pid, acc, sites, useq)
            except Exception:
                continue
            if r and r.get("usable"):
                cands.append(r)
                print(f"    ✓ {pid}  链{r['chain']} cov={r['cov']:.0%} "
                      f"范围{r['range']} res={r['resolution']} 覆盖位点{r['sites_covered']}")
        cands.sort(key=lambda x: (x["resolution"] if x["resolution"] else 99))
        report[name] = {"uniprot": acc, "sites": sites,
                        "n_entries": total, "candidates": cands}
        if not cands:
            print("    ✗ 未找到符合条件的非磷酸化结构")

    outp = OUT / "unphos_candidates.json"
    merged = {}
    if outp.exists():
        try: merged = json.loads(outp.read_text())
        except Exception: pass
    merged.update(report)
    outp.write_text(json.dumps(merged, indent=2, ensure_ascii=False))
    print(f"\n[OK] {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
