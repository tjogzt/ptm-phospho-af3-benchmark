# 第二批作业上传指南（prkaca / egfr / lyz / polr2a）

生成时间：2026-10-23 · 来源：`audit/af3_tosubmit.json` · 校验：`audit/site_recheck_tosubmit.md`

## 一、文件

| 文件 | 内容 |
|---|---|
| `af3_jobs/batches/batch_tosubmit_1.json` | prkaca + egfr，共 12 个作业 |
| `af3_jobs/batches/batch_tosubmit_2.json` | lyz + polr2a，共 12 个作业 |
| `af3_jobs/tosubmit_v3_jobs.json` | 全部 24 个作业（如需一次提交） |

**为什么是 24 个而不是 8 个**：`scripts/af3/download_af3_results.py` 用正则 `^(\w+)_(phos|unphos)_seed(\d+)$` 解析作业名，**一个作业只能含一个 seed**。原 `af3_tosubmit.json` 里每个作业带 3 个 seed，直接上传会导致解析失败。已按 4 蛋白 × 2 状态 × 3 seed 拆开，与历史已成功的 88 个作业完全一致。

## 二、格式已逐字段比对本批次与历史成功批次

| 项 | 本批 | 历史批次 | 一致 |
|---|---|---|---|
| 顶层结构 | 数组 `[ {...} ]` | 同 | ✓ |
| job 顶层字段 | `name / modelSeeds / sequences / dialect / version` | 同 | ✓ |
| `dialect` / `version` | `alphafoldserver` / `3` | 同 | ✓ |
| `modelSeeds` | `["1"]`（字符串数组） | 同 | ✓ |
| proteinChain 字段 | `sequence / count / useStructureTemplate` | 同 | ✓ |
| 修饰 | `{"ptmType":"CCD_*","ptmPosition":N}` | 同 | ✓ |

## 三、上传步骤

1. 打开 https://alphafoldserver.com ，用 Google 账号登录。
2. 选「Upload JSON」→ 先传 `batch_tosubmit_1.json`（12 个）。
3. 上传后核对页面显示的修饰：
   - `prkaca_phos_seed1`：CCD_TPO @ **198**
   - `egfr_phos_seed1`：CCD_PTR @ **869**
   - `lyz_phos_seed1`：CCD_PTR @ **41**
   - `polr2a_phos_seed1`：CCD_SEP @ **1843** + **1845**
4. 提交，等状态 `Running → Complete`。
5. 次日再传 `batch_tosubmit_2.json`（配额通常 20–30 job/天，按 UTC 重置）。

## 四、提交前必读：两个冲突

1. **POLR2A 位点已改**：旧批次用的是 1614/1621，本批改为 1843/1845（UniProt 实测两者均为 Phosphoserine）。作业名相同 → 旧作业结果**不要再下载**，否则与本次同名冲突。
2. **PRKACA 可能已提交过**：旧批次 prkaca 的序列(351)与位点(198)与本批完全相同，但 `af3_ledger.json` 里没有 prkaca 条目（当时结果未下载）。建议先去网页看 7 月那批 prkaca 是否还能下载；能下载就不必重传，可省 6 个配额。若查不到，直接传即可（内容一致，无风险）。

## 五、下载与后续

1. 下载 zip 放入 `results/af3/downloads/`（文件名保持 `folds_*.zip`）。
2. 我这边执行：

```bash
python3 scripts/af3/download_af3_results.py --extract      # 幂等解析
python3 scripts/analysis/analyze_af3_ptm_effect.py          # 重跑 ΔpLDDT / ΔpTM
python3 scripts/analysis/compare_af3_vs_pdb.py --all-models # 重跑 AF3 vs 实验位移
python3 scripts/analysis/plot_af3_vs_pdb.py
```

3. 解析后请抽查：`af3_summary.csv` 中 polr2a 的修饰位点应为 **1843/1845**，若为 1614/1621 说明混入旧结果，需剔除。
