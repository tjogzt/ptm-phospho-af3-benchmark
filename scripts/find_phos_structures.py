#!/usr/bin/env python3
"""
find_phos_structures.py — 按"真实磷酸化修饰"检索候选 PDB 结构
================================================================
动机：原 benchmark 按蛋白名检索 PDB 再人工填 ID，导致三类错误混入
（PDB ID 张冠李戴、非磷酸化结构被当磷酸化态、结构域错配）。
本脚本改为以**修饰证据**为准：检索该 UniProt 的全部 PDB entry，
逐一取文件头部的 MODRES/SEQADV 记录，只保留真正含 SEP/TPO/PTR 者。

实现要点：
  - RCSB Search API 按 UniProt accession 列出 entry
  - 用 HTTP Range 只下载每份 PDB 的前 256 KB（MODRES/SEQRES 均在头部），
    避免整库下载
  - 解析 MODRES 得到"哪条链、哪个残基、何种修饰"
  - 额外报告分辨率与标题，便于人工裁决

输出：audit/phos_structure_candidates.json + 控制台表格

用法：
  python3 scripts/find_phos_structures.py                # 全部待修复蛋白
  python3 scripts/find_phos_structures.py --only CDK7 IRF3
"""
import sys, json, argparse, time
from pathlib import Path
import urllib.request

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

OUT = _ROOT / "audit"; OUT.mkdir(exist_ok=True)
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
FILES = "https://files.rcsb.org/download/{pid}.pdb"
PHOS = {"SEP", "TPO", "PTR", "NEP"}
HEADER_BYTES = 262144


def search_entries(acc, rows=200):
    q = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers"
                         ".reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows},
                            "results_content_type": ["experimental"]},
    }
    req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(),
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=40) as r:
            d = json.loads(r.read().decode())
        return [x["identifier"] for x in d.get("result_set", [])], d.get("total_count", 0)
    except Exception as e:
        print(f"    [SEARCH FAIL] {acc}: {type(e).__name__}")
        return [], 0


def fetch_header(pid):
    """Range 请求取 PDB 头部（含 MODRES / TITLE / REMARK 2 分辨率）"""
    req = urllib.request.Request(FILES.format(pid=pid),
                                 headers={"Range": f"bytes=0-{HEADER_BYTES}"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.read().decode("utf-8", "ignore")
    except Exception:
        return ""


def parse_header(txt):
    """返回 (title, modres[list], resolution)"""
    title, mods, res = [], [], None
    for line in txt.splitlines():
        rec = line[:6].strip()
        if rec == "TITLE":
            title.append(line[10:].strip())
        elif rec == "MODRES":
            resn = line[12:15].strip()
            if resn in PHOS:
                mods.append({"resname": resn, "chain": line[16],
                             "seq": int(line[18:22].strip() or 0),
                             "std": line[24:27].strip()})
        elif rec == "REMARK" and line[7:10].strip() == "2" and "RESOLUTION" in line.upper():
            if res is None:
                parts = line[10:].replace("ANGSTROM.", "").replace("ANGSTROM", "").split()
                for p in parts:
                    try:
                        v = float(p)
                        if 0.5 <= v <= 15:
                            res = v; break
                    except ValueError:
                        continue
    return " ".join(title), mods, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None, help="只查指定蛋白（benchmark key）")
    ap.add_argument("--delay", type=float, default=0.15)
    args = ap.parse_args()

    # 需要重新检索的蛋白：phos 结构经审计不成立者 + 单结构组
    TARGETS = ["CDK7", "IRF3", "P38alpha", "Hsp27", "c-Src_kinase", "CFTR_R_domain",
               "EGFR_kinase", "AMPK_alpha2", "PhK_gamma", "PKA_Calpha"]
    targets = args.only if args.only else TARGETS

    report = {}
    print("="*100)
    print("  按真实磷酸化修饰检索候选结构")
    print("="*100)

    for name in targets:
        info = BENCHMARK_PROTEINS.get(name)
        if not info:
            print(f"[SKIP] 未知蛋白 {name}"); continue
        acc = info["uniprot"]; gene = info["gene"]
        sites = info.get("phos_sites", [])
        print(f"\n### {name} ({gene} / {acc}) 目标位点: {sites}")
        entries, total = search_entries(acc)
        print(f"    UniProt 关联 entry: {total}（取前 {len(entries)} 个验证）")

        hits = []
        for i, pid in enumerate(entries, 1):
            txt = fetch_header(pid)
            if not txt:
                continue
            title, mods, res = parse_header(txt)
            if not mods:
                continue
            # 只接受链匹配或任意链的磷酸化（记录链信息供人工判断）
            hit = {"pdb": pid, "title": title[:130], "resolution": res,
                   "modifications": mods,
                   "n_phos": len(mods),
                   "chains": sorted({m["chain"] for m in mods})}
            hits.append(hit)
            print(f"    ✓ {pid}  res={res}  n_phos={len(mods)}  chains={hit['chains']}  "
                  f"seqs={sorted({m['seq'] for m in mods})[:8]}")
            time.sleep(args.delay)

        # 按磷酸化位点数与分辨率排序
        hits.sort(key=lambda h: (-h["n_phos"], h["resolution"] if h["resolution"] else 99))
        report[name] = {"gene": gene, "uniprot": acc, "declared_sites": sites,
                        "n_entries": total, "candidates": hits}
        if not hits:
            print("    ✗ 未找到含 SEP/TPO/PTR 的实验结构")

    out_p = OUT/"phos_structure_candidates.json"
    merged = {}
    if out_p.exists():
        try: merged = json.loads(out_p.read_text())
        except Exception: merged = {}
    merged.update(report)                    # 增量合并
    out_p.write_text(json.dumps(merged, indent=2, ensure_ascii=False))
    report = merged

    print("\n" + "="*100)
    print("  汇总：各蛋白真实磷酸化候选结构数")
    print("="*100)
    for name, r in report.items():
        c = r["candidates"]
        top = c[0]["pdb"] if c else "—"
        print(f"  {name:<22} 候选 {len(c):>3} 个   最佳: {top}")
    print(f"\n[OK] {OUT}/phos_structure_candidates.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
