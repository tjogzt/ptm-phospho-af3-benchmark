# af3_tosubmit.json 位点复核（提交前）

- 输入：`audit/af3_tosubmit.json`（8 个作业组）
- 参照：`data/curated/uniprot_sequences.json`（本地 UniProt canonical 序列）
- 拆分产物：`af3_jobs/tosubmit_v3_jobs.json`（共 24 个作业，2 批）

## 逐作业校验

| 作业组 | UniProt | 长度(声明/实际) | 与UniProt缓存 | 位点(UniProt) | 期望残基 | 实际残基 | OK |
|---|---|---|---|---|---|---|---|
| prkaca_unphos | P17612 | 351/351 | 一致 | 198:T | T | ✓ | ✓ |
| prkaca_phos | P17612 | 351/351 | 一致 | 198:T | T | ✓ | ✓ |
| egfr_unphos | P00533 | 1210/1210 | 一致 | 869:Y | Y | ✓ | ✓ |
| egfr_phos | P00533 | 1210/1210 | 一致 | 869:Y | Y | ✓ | ✓ |
| lyz_unphos | P00698 | 147/147 | 一致 | 41:Y | Y | ✓ | ✓ |
| lyz_phos | P00698 | 147/147 | 一致 | 41:Y | Y | ✓ | ✓ |
| polr2a_unphos | P24928 | 1970/1970 | 一致 | 1843:S, 1845:S | S, S | ✓ | ✓ |
| polr2a_phos | P24928 | 1970/1970 | 一致 | 1843:S, 1845:S | S, S | ✓ | ✓ |

## 位点裁决依据（原文摘录）

- **prkaca_unphos**（P17612）：benchmark Thr197 是成熟编号；UniProt[197]=W，真实 p-site 为 196/198，取 198（激活环主位点）
- **prkaca_phos**（P17612）：benchmark Thr197 是成熟编号；UniProt[197]=W，真实 p-site 为 196/198，取 198（激活环主位点）
- **egfr_unphos**（P00533）：benchmark Tyr845 是成熟编号(+24 信号肽)；UniProt[869] 为注释 pTyr
- **egfr_phos**（P00533）：benchmark Tyr845 是成熟编号(+24 信号肽)；UniProt[869] 为注释 pTyr
- **lyz_unphos**（P00698）：人工引入位点；成熟 Tyr23 + 18(信号肽) = UniProt[41]=Y；UniProt 无天然磷酸化注释（阴性对照符合设计）
- **lyz_phos**（P00698）：人工引入位点；成熟 Tyr23 + 18(信号肽) = UniProt[41]=Y；UniProt 无天然磷酸化注释（阴性对照符合设计）
- **polr2a_unphos**（P24928）：benchmark Ser2/Ser5 是 CTD heptad 内相对编号；映射到 CTD 区已注释 pSer 1843/1845（YSPTSPS 第2/第5位）
- **polr2a_phos**（P24928）：benchmark Ser2/Ser5 是 CTD heptad 内相对编号；映射到 CTD 区已注释 pSer 1843/1845（YSPTSPS 第2/第5位）

## 问题清单

- 无。长度、序列、位点残基、修饰类型四项全部通过。

## UniProt 权威注释复核（2026-10-23 实查 REST）

来源：`https://rest.uniprot.org/uniprotkb/{ACC}.json?fields=ft_mod_res`，原始证据存 `audit/uniprot_site_evidence.json`

| 蛋白 | UniProt | 提交位点 | UniProt 注释 | 结论 |
|---|---|---|---|---|
| PRKACA | P17612 | T198 | `Phosphothreonine; by PDPK1` | ✓ 激活环主位点（PDK1）。注意 196 位亦注释为 Phosphothreonine（次级位点） |
| EGFR | P00533 | Y869 | `Phosphotyrosine` | ✓ |
| LYZ | P00698 | Y41 | 该蛋白 **0 条** modified residue 注释 | ✓ 确为人工引入位点，符合阴性对照设计 |
| POLR2A | P24928 | S1843 / S1845 | 两者均为 `Phosphoserine` | ✓ CTD heptad 区；邻位 1840/1847/1849/1850 亦有注释 |

## 与历史提交的冲突（提交前必须处理）

| 问题 | 详情 | 建议 |
|---|---|---|
| POLR2A 位点已改 | `af3_jobs/af3_jobs_v2.json` 中 polr2a 用的是 **1614/1621**，本次改为 **1843/1845**；但作业名同为 `polr2a_phos_seed{1,2,3}` | 旧作业结果**切勿再下载**（若日后下载会与本批同名冲突）。新结果入账时核对 `modifications` 位点应为 1843/1845 |
| PRKACA 重复提交风险 | v2 批次中 prkaca 已提交过，**序列(351)与位点(198)与本次完全一致**；`af3_ledger.json` 中无 prkaca 条目 → 说明当时结果未下载 | 先在 AlphaFold Server 网页查 7 月那批 prkaca 是否还能下载；能下载就不必重传，省 6 个配额 |

## 批次

- `af3_jobs/batches/batch_tosubmit_1.json`：12 个作业
- `af3_jobs/batches/batch_tosubmit_2.json`：12 个作业
