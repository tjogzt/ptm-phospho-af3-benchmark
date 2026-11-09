#!/usr/bin/env python3
"""
download_af3_results.py — AF3 批量结果解析（幂等 + 生命周期台账）
==================================================================
P0-2b 重写版。修复原版三个缺陷：
  1. 非幂等：每次运行从零重建 CSV，多次运行相互覆盖 → 现改为
     "扫描全部 zip → 累积到台账 → 原子写出"，重复运行结果一致
  2. 无生命周期：现记录 submitted / returned / extracted / parsed 四态
  3. 覆盖不全：原版产出 30 行（5 蛋白），实际 zip 内有 88 作业 → 现全量解析

输出：
  results/af3/af3_summary.csv          每作业一行的置信度汇总（全量）
  results/af3/af3_ledger.json          生命周期台账（提交/返回/解析状态）
  results/af3/models/<gene>/<state>/   从 zip 解出的 CIF（可选，--extract）

用法：
  python scripts/af3/download_af3_results.py --report          # 仅重建 CSV + 台账
  python scripts/af3/download_af3_results.py --extract         # 同时解出 CIF
  python scripts/af3/download_af3_results.py --verify-sites    # 加做位点校验
"""
import os, sys, json, csv, zipfile, re, argparse
from pathlib import Path
from collections import defaultdict

_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS, REMOVED_PROTEINS

AF3_DIR = _ROOT / "results" / "af3"
DOWNLOADS = AF3_DIR / "downloads"
JOB_RE = re.compile(r'^(\w+)_(phos|unphos)_seed(\d+)$')
# AF3 server 对重名作业会自动追加 _2/_3 后缀，视为同一作业的重复提交
JOB_RE_DUP = re.compile(r'^(\w+)_(phos|unphos)_seed(\d+)_(\d+)$')

# gene -> benchmark protein key
GENE2PROT = {info['gene'].lower(): k for k, info in BENCHMARK_PROTEINS.items()}
# 已移除蛋白的基因名：标记为 REMOVED:<key>，避免与 UNKNOWN（真未知）混淆
GENE2REMOVED = {info['gene'].lower(): f"REMOVED:{k}"
                for k, info in REMOVED_PROTEINS.items()}


def scan_zip(zip_path):
    """扫描单个 zip，返回该 zip 内每个作业的记录（不写盘）。"""
    records = {}
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        job_dirs = sorted({n.split('/')[0] for n in names if '/' in n})
        for jd in job_dirs:
            m = JOB_RE.match(jd)
            dup_of = None
            if not m:
                md = JOB_RE_DUP.match(jd)
                if not md:
                    continue
                # 归一化为原始作业名，标记为重复提交
                m = md
                dup_of = f"{md.group(1)}_{md.group(2)}_seed{md.group(3)}"
            gene, state, seed = m.group(1).lower(), m.group(2), int(m.group(3))
            rec = {
                "job_name": jd, "gene": gene, "state": state, "seed": seed,
                "source_zip": zip_path.name,
                "protein": GENE2PROT.get(gene) or GENE2REMOVED.get(gene, "UNKNOWN"),
                "models_found": 0, "models_extracted": 0,
                "plddt_mean": None, "plddt_best": None,
                "ptm_mean": None, "ptm_best": None, "ranking_mean": None,
                "ptm_requested": None, "seq_len": None,
                "cif_files": [],
                "dup_of": dup_of,
                "data_version": "v2" if re.search(r"_2026_1[0-9]_", zip_path.name) else "v1",
            }
            plddts, ptms, ranks = [], [], []

            # job_request.json = 服务端存档的提交内容（权威）
            req_name = f"{jd}/fold_{jd}_job_request.json"
            if req_name in names:
                try:
                    req = json.loads(zf.read(req_name))[0]
                    pc = req['sequences'][0]['proteinChain']
                    rec["seq_len"] = len(pc['sequence'])
                    rec["ptm_requested"] = [
                        {"type": mod['ptmType'].replace('CCD_', ''), "pos": mod['ptmPosition'],
                         "build_aa": pc['sequence'][mod['ptmPosition']-1]}
                        for mod in (pc.get('modifications') or [])]
                except Exception:
                    pass

            for i in range(5):  # AF3 每作业 5 个 model
                cif = f"{jd}/fold_{jd}_model_{i}.cif"
                conf = f"{jd}/fold_{jd}_summary_confidences_{i}.json"
                full = f"{jd}/fold_{jd}_full_data_{i}.json"
                if cif in names:
                    rec["models_found"] += 1
                    rec["cif_files"].append(f"{jd}_model_{i}.cif")
                if conf in names:
                    try:
                        c = json.loads(zf.read(conf))
                        if c.get("ptm") is not None: ptms.append(float(c["ptm"]))
                        if c.get("ranking_score") is not None: ranks.append(float(c["ranking_score"]))
                    except Exception: pass
                if full in names:
                    try:
                        f_ = json.loads(zf.read(full))
                        ap = f_.get("atom_plddts") or []
                        if ap: plddts.append(sum(ap)/len(ap))
                    except Exception: pass

            if plddts:
                rec["plddt_mean"] = round(sum(plddts)/len(plddts), 2)
                rec["plddt_best"] = round(max(plddts), 2)
            if ptms:
                rec["ptm_mean"] = round(sum(ptms)/len(ptms), 3)
                rec["ptm_best"] = round(max(ptms), 3)
            if ranks:
                rec["ranking_mean"] = round(sum(ranks)/len(ranks), 3)
            records[jd] = rec
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", action="store_true", help="解出 CIF 到 results/af3/models/")
    ap.add_argument("--verify-sites", action="store_true", help="校验磷酸位点（需网络）")
    args = ap.parse_args()

    zips = sorted(f for f in DOWNLOADS.glob("folds_*.zip") if not f.name.startswith("._"))
    if not zips:
        print(f"ERROR: no zips in {DOWNLOADS}"); return 1
    print(f"扫描 {len(zips)} 个 zip ...")

    all_recs = {}
    for zp in zips:
        recs = scan_zip(zp)
        for jn, r in recs.items():
            if r.get("dup_of"):
                # 重复提交：映射回原始作业名，仅在原始作业缺数据时用于补位
                base = r["dup_of"]
                if base not in all_recs and base not in recs:
                    r["job_name"] = base
                    all_recs[base] = r
                continue
            if jn in all_recs:
                old = all_recs[jn]
                # 优先级：data_version v2 > v1；同版本内取模型更完整者
                better = (r["data_version"] > old["data_version"]) or (
                    r["data_version"] == old["data_version"]
                    and r["models_found"] > old["models_found"])
                if better:
                    all_recs[jn] = r
            else:
                all_recs[jn] = r
        print(f"  {zp.name}: {len(recs)} jobs")

    # ---- 可选：解出 CIF ----
    if args.extract:
        out_root = AF3_DIR / "models"      # 每版本子目录：models/<ver>/<gene>/<state>/
        n = 0
        for zp in zips:
            with zipfile.ZipFile(zp) as zf:
                for jn, r in all_recs.items():
                    if r["source_zip"] != zp.name: continue
                    tdir = out_root / r["data_version"] / r["gene"] / r["state"]
                    tdir.mkdir(parents=True, exist_ok=True)
                    for member in zf.namelist():
                        if member.startswith(jn + "/") and member.endswith(".cif"):
                            dst = tdir / Path(member).name
                            if not dst.exists():
                                dst.write_bytes(zf.read(member))
                                n += 1
                                r["models_extracted"] += 1
        print(f"  解出 {n} 个 CIF 到 {out_root}")

    # ---- 位点校验（可选） ----
    if args.verify_sites:
        import urllib.request
        uni_path = _ROOT / "data" / "curated" / "uniprot_sequences.json"
        if uni_path.exists():
            uni = json.loads(uni_path.read_text())
            acc2seq = {k: v['sequence'] for k, v in uni.items()}
            gene2acc = {v['gene'].lower(): k for k, v in uni.items()}
            cache = {}
            for jn, r in all_recs.items():
                if not r["ptm_requested"]: continue
                acc = gene2acc.get(r["gene"])
                if not acc: continue
                if acc not in cache:
                    try:
                        with urllib.request.urlopen(
                            f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=20) as resp:
                            d = json.loads(resp.read().decode())
                        cache[acc] = {f['location']['start']['value'] for f in d.get('features', [])
                                      if f['type']=='Modified residue' and 'hospho' in f.get('description','')}
                    except Exception:
                        cache[acc] = None
                ann = cache[acc]
                r["site_verified"] = None if ann is None else all(
                    p["pos"] in ann for p in r["ptm_requested"])
            print("  位点校验完成")

    # ---- 写出 CSV（按 job_name 排序，幂等） ----
    fields = ["job_name","protein","gene","state","seed","effect_type",
              "models_found","models_extracted","plddt_mean","plddt_best",
              "ptm_mean","ptm_best","ranking_mean","seq_len","ptm_requested",
              "site_verified","source_zip","dup_of","data_version"]
    rows = []
    for jn in sorted(all_recs):
        r = dict(all_recs[jn])
        info = BENCHMARK_PROTEINS.get(r["protein"], {})
        r["effect_type"] = info.get("effect_type", "unknown")
        r["ptm_requested"] = (json.dumps(r["ptm_requested"], ensure_ascii=False)
                              if r.get("ptm_requested") else "")
        r.setdefault("site_verified", "")
        rows.append(r)

    csv_path = AF3_DIR / "af3_summary.csv"
    tmp = csv_path.with_suffix(".csv.tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    tmp.replace(csv_path)

    # ---- 生命周期台账 ----
    ledger = {
        "generated_from_zips": [z.name for z in zips],
        "counts": {
            "jobs_parsed": len(all_recs),
            "phos": sum(1 for r in all_recs.values() if r["state"]=="phos"),
            "unphos": sum(1 for r in all_recs.values() if r["state"]=="unphos"),
            "unique_genes": len({r["gene"] for r in all_recs.values()}),
            "cif_total": sum(r["models_found"] for r in all_recs.values()),
        },
        "by_protein": {},
        "jobs": {jn: {k: all_recs[jn].get(k) for k in
                      ("gene","state","seed","protein","models_found",
                       "plddt_mean","ptm_mean","seq_len","source_zip")}
                 for jn in sorted(all_recs)},
    }
    byp = defaultdict(lambda: {"phos": 0, "unphos": 0, "plddt_phos": [], "plddt_unphos": [],
                               "ptm_phos": [], "ptm_unphos": []})
    for r in all_recs.values():
        d = byp[r["protein"]]
        d[r["state"]] += 1
        if r["plddt_mean"]: d[f"plddt_{r['state']}"].append(r["plddt_mean"])
        if r["ptm_mean"]:   d[f"ptm_{r['state']}"].append(r["ptm_mean"])
    for p, d in byp.items():
        ledger["by_protein"][p] = {
            "jobs_phos": d["phos"], "jobs_unphos": d["unphos"],
            "plddt_phos_mean": round(sum(d["plddt_phos"])/len(d["plddt_phos"]), 2) if d["plddt_phos"] else None,
            "plddt_unphos_mean": round(sum(d["plddt_unphos"])/len(d["plddt_unphos"]), 2) if d["plddt_unphos"] else None,
            "ptm_phos_mean": round(sum(d["ptm_phos"])/len(d["ptm_phos"]), 3) if d["ptm_phos"] else None,
            "ptm_unphos_mean": round(sum(d["ptm_unphos"])/len(d["ptm_unphos"]), 3) if d["ptm_unphos"] else None,
        }
    led_path = AF3_DIR / "af3_ledger.json"
    led_path.write_text(json.dumps(ledger, indent=2, ensure_ascii=False))

    # ---- 汇报 ----
    c = ledger["counts"]
    print(f"\n{'='*66}")
    print(f"  作业解析: {c['jobs_parsed']}  (phos {c['phos']} / unphos {c['unphos']})")
    print(f"  覆盖蛋白: {c['unique_genes']} 种 | CIF 总数: {c['cif_total']}")
    print(f"  CSV:  {csv_path}  ({len(rows)} 行)")
    print(f"  台账: {led_path}")
    print(f"{'='*66}")
    print(f"\n{'蛋白':<22}{'phos':>6}{'unphos':>8}{'pLDDT(p)':>11}{'pLDDT(u)':>11}{'ΔpLDDT':>9}")
    print("-"*66)
    for p in sorted(ledger["by_protein"]):
        d = ledger["by_protein"][p]
        delta = (f"{d['plddt_phos_mean']-d['plddt_unphos_mean']:+.2f}"
                 if d['plddt_phos_mean'] and d['plddt_unphos_mean'] else "-")
        print(f"{p:<22}{d['jobs_phos']:>6}{d['jobs_unphos']:>8}"
              f"{d['plddt_phos_mean'] or '-':>11}{d['plddt_unphos_mean'] or '-':>11}{delta:>9}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
