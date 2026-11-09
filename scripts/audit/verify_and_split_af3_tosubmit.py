#!/usr/bin/env python3
"""校验 af3_tosubmit.json 的位点映射，并拆分为 AlphaFold Server 可上传的作业数组。

用法:
    python3 scripts/audit/verify_and_split_af3_tosubmit.py

输入:  audit/af3_tosubmit.json            （位点已裁决的待提交清单）
       data/curated/uniprot_sequences.json（本地 UniProt 权威序列缓存）
输出:  af3_jobs/tosubmit_v3_jobs.json     （全部作业，数组，dialect=alphafoldserver/version=3）
       af3_jobs/batches/batch_tosubmit_1.json / _2.json （按每日配额分批）
       audit/site_recheck_tosubmit.json   （逐作业校验明细）
       audit/site_recheck_tosubmit.md     （人读校验报告）

命名规范必须与 scripts/af3/download_af3_results.py 的正则一致:
    JOB_RE = ^(\\w+)_(phos|unphos)_seed(\\d+)$   →  每个作业只含 1 个 seed
"""

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / "audit" / "af3_tosubmit.json"
UNIPROT = ROOT / "data" / "curated" / "uniprot_sequences.json"
OUT_JOBS = ROOT / "af3_jobs" / "tosubmit_v3_jobs.json"
OUT_BATCH_DIR = ROOT / "af3_jobs" / "batches"
OUT_JSON = ROOT / "audit" / "site_recheck_tosubmit.json"
OUT_MD = ROOT / "audit" / "site_recheck_tosubmit.md"

BATCH_SIZE = 12  # AlphaFold Server 每日配额通常 20-30 job/天，保守取 12


def die(msg):
    print(f"[FATAL] {msg}", file=sys.stderr)
    sys.exit(1)


def main():
    for p in (SRC, UNIPROT):
        if not p.exists():
            die(f"缺少输入文件: {p}")

    tosubmit = json.loads(SRC.read_text(encoding="utf-8"))
    uniprot = json.loads(UNIPROT.read_text(encoding="utf-8"))

    jobs_src = tosubmit["jobs"]
    print(f"待提交作业组: {len(jobs_src)}（来自 {SRC.name}）")

    report = []
    out_jobs = []
    problems = []

    for name, rec in jobs_src.items():
        job = rec["job"]
        acc = rec["uniprot"]
        seq = job["sequences"][0]["proteinChain"]["sequence"]
        mods = job["sequences"][0]["proteinChain"].get("modifications", [])
        sites = rec["sites_uniprot"]
        expect_res = rec["site_residues"]
        gene = name.split("_")[0]  # prkaca / egfr / lyz / polr2a

        # --- 校验 1: 长度 ---
        len_ok = len(seq) == rec["seq_len"]
        if not len_ok:
            problems.append(f"{name}: 序列长度 {len(seq)} != 声明 {rec['seq_len']}")

        # --- 校验 2: 与本地 UniProt 缓存逐字符比对 ---
        ref = uniprot.get(acc, {}).get("sequence")
        if ref is None:
            seq_match, diff_note = None, f"本地缓存无 {acc}，无法比对"
            problems.append(f"{name}: 本地 uniprot_sequences.json 缺少 {acc}")
        elif ref == seq:
            seq_match, diff_note = True, "与本地 UniProt 缓存完全一致"
        else:
            seq_match = False
            diffs = [
                (i + 1, a, b)
                for i, (a, b) in enumerate(zip(ref, seq))
                if a != b
            ]
            if len(ref) != len(seq):
                diff_note = f"长度不同 ref={len(ref)} vs job={len(seq)}"
            else:
                diff_note = f"{len(diffs)} 处残基差异，前 5 处: " + "; ".join(
                    f"{i}:{a}->{b}" for i, a, b in diffs[:5]
                )
            problems.append(f"{name}: 与 UniProt 缓存序列不一致 → {diff_note}")

        # --- 校验 3: 位点残基 ---
        site_checks = []
        for pos, exp in zip(sites, expect_res):
            actual = seq[pos - 1] if pos - 1 < len(seq) else "?"
            ok = actual == exp
            site_checks.append(
                {"pos_uniprot": pos, "expect": exp, "actual": actual, "ok": ok}
            )
            if not ok:
                problems.append(f"{name}: 位点 {pos} 期望 {exp}，实际 {actual}")

        # --- 校验 4: 修饰条目与位点一致（unphos 应无修饰，phos 应覆盖全部位点）---
        state = "phos" if name.endswith("_phos") else "unphos"
        mod_positions = sorted(m["ptmPosition"] for m in mods)
        if state == "unphos":
            mod_ok = not mods
            if not mod_ok:
                problems.append(f"{name}: unphos 作业不应带修饰，实际 {mod_positions}")
        else:
            mod_ok = mod_positions == sorted(sites)
            if not mod_ok:
                problems.append(
                    f"{name}: modifications 位点 {mod_positions} != sites_uniprot {sorted(sites)}"
                )

        # --- 修饰类型与残基是否匹配（CCD 码 vs 残基） ---
        ccd_expect = {"T": "CCD_TPO", "Y": "CCD_PTR", "S": "CCD_SEP"}
        for m in mods:
            pos = m["ptmPosition"]
            res = seq[pos - 1] if pos - 1 < len(seq) else "?"
            if ccd_expect.get(res) != m["ptmType"]:
                problems.append(
                    f"{name}: 位点 {pos} 残基 {res} 与 {m['ptmType']} 不匹配，应为 {ccd_expect.get(res)}"
                )

        # --- 拆分为每 seed 一个作业（与 download_af3_results.py 正则一致） ---
        for seed in job["modelSeeds"]:
            chain = {
                "sequence": seq,
                "count": 1,
                "useStructureTemplate": True,
            }
            if mods:
                chain["modifications"] = mods
            out_jobs.append(
                {
                    "name": f"{gene}_{state}_seed{seed}",
                    "modelSeeds": [str(seed)],
                    "sequences": [{"proteinChain": chain}],
                    "dialect": "alphafoldserver",
                    "version": 3,
                }
            )

        report.append(
            {
                "job_group": name,
                "uniprot": acc,
                "gene": gene,
                "state": state,
                "seq_len_declared": rec["seq_len"],
                "seq_len_actual": len(seq),
                "len_ok": len_ok,
                "identical_to_uniprot_cache": seq_match,
                "seq_diff_note": diff_note,
                "sites": site_checks,
                "modifications": mods,
                "mod_positions": mod_positions,
                "mod_matches_sites": mod_ok,
                "mapping_note": rec.get("mapping_note", ""),
                "effect_type": rec.get("effect_type", ""),
            }
        )

    # --- 写产物 ---
    OUT_JOBS.parent.mkdir(parents=True, exist_ok=True)
    OUT_JOBS.write_text(
        json.dumps(out_jobs, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    OUT_BATCH_DIR.mkdir(parents=True, exist_ok=True)
    batches = [
        out_jobs[i : i + BATCH_SIZE] for i in range(0, len(out_jobs), BATCH_SIZE)
    ]
    batch_index = []
    for i, b in enumerate(batches, 1):
        p = OUT_BATCH_DIR / f"batch_tosubmit_{i}.json"
        p.write_text(json.dumps(b, ensure_ascii=False, indent=2), encoding="utf-8")
        batch_index.append(
            {"batch": i, "file": str(p.relative_to(ROOT)), "n_jobs": len(b)}
        )

    OUT_JSON.write_text(
        json.dumps(
            {"n_job_groups": len(report), "checks": report, "problems": problems},
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    # --- Markdown 报告 ---
    lines = [
        "# af3_tosubmit.json 位点复核（提交前）",
        "",
        f"- 输入：`audit/af3_tosubmit.json`（{len(report)} 个作业组）",
        "- 参照：`data/curated/uniprot_sequences.json`（本地 UniProt canonical 序列）",
        f"- 拆分产物：`af3_jobs/tosubmit_v3_jobs.json`（共 {len(out_jobs)} 个作业，{len(batches)} 批）",
        "",
        "## 逐作业校验",
        "",
        "| 作业组 | UniProt | 长度(声明/实际) | 与UniProt缓存 | 位点(UniProt) | 期望残基 | 实际残基 | OK |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in report:
        sites_s = ", ".join(
            f"{c['pos_uniprot']}:{c['actual']}"
            + ("" if c["ok"] else f"(期望{c['expect']})")
            for c in r["sites"]
        )
        exp_s = ", ".join(str(c["expect"]) for c in r["sites"])
        ok_s = "✓" if all(c["ok"] for c in r["sites"]) and r["len_ok"] else "✗"
        lines.append(
            f"| {r['job_group']} | {r['uniprot']} | {r['seq_len_declared']}/{r['seq_len_actual']} "
            f"| {'一致' if r['identical_to_uniprot_cache'] else r['seq_diff_note']} "
            f"| {sites_s} | {exp_s} | {'✓' if all(c['ok'] for c in r['sites']) else '✗'} | {ok_s} |"
        )
    lines += ["", "## 位点裁决依据（原文摘录）", ""]
    for r in report:
        lines.append(f"- **{r['job_group']}**（{r['uniprot']}）：{r['mapping_note']}")
    lines += ["", "## 问题清单", ""]
    if problems:
        for p in problems:
            lines.append(f"- ⚠️ {p}")
    else:
        lines.append("- 无。长度、序列、位点残基、修饰类型四项全部通过。")
    lines += ["", "## 批次", ""]
    for b in batch_index:
        lines.append(f"- `{b['file']}`：{b['n_jobs']} 个作业")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # --- 控制台汇总 ---
    print("\n=== 校验结果 ===")
    for r in report:
        flag = "OK " if all(c["ok"] for c in r["sites"]) and r["len_ok"] else "BAD"
        print(
            f"[{flag}] {r['job_group']:16s} {r['uniprot']} len={r['seq_len_actual']} "
            f"uniprot={'一致' if r['identical_to_uniprot_cache'] else '不一致'} "
            + " ".join(
                f"{c['pos_uniprot']}{c['actual']}{'✓' if c['ok'] else '✗'}"
                for c in r["sites"]
            )
        )
    print(f"\n拆分为 {len(out_jobs)} 个作业 → {OUT_JOBS.relative_to(ROOT)}")
    for b in batch_index:
        print(f"  批次 {b['batch']}: {b['file']} ({b['n_jobs']} jobs)")
    if problems:
        print(f"\n⚠️ {len(problems)} 处问题：")
        for p in problems:
            print(f"   - {p}")
        sys.exit(1)
    print("\n全部通过。")


if __name__ == "__main__":
    main()
