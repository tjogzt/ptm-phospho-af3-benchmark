# AlphaFold Server 上传操作指南

## 一、准备

作业文件已就绪，无需再生成：

```
af3_jobs/
├── batched/batch_1.json   ← 第 1 批（20 作业）
├── batched/batch_2.json   ← 第 2 批（20 作业）
├── batched/batch_3.json   ← 第 3 批（20 作业）
├── batched/batch_index.json  批次清单
├── af3_jobs_v2.json       ← 全部 60 作业（如需一次性提交）
└── manifest_v2.json       ← 蛋白/位点对照表
```

**分 3 批的原因**：AlphaFold Server 对每个账号有**每日作业配额**（通常 20–30 个 job/天）。
一次上传 60 个会被部分拒绝，且失败原因不易定位。

**批次顺序已按优先级排好**：成对蛋白（能做 AF3-vs-实验验证的那 6 组）的 phos 作业排在最前。

## 二、上传步骤

1. **打开** https://alphafoldserver.com
2. **登录**：用 Google 账号（必须是 Google 账号，不支持邮箱密码）
3. **接受条款**：首次使用需同意 AlphaFold Server 使用条款
4. **上传 JSON**：
   - 页面顶部/左侧有「Upload」或「+ New job」→ 选择「Upload JSON」
   - 选中 `af3_jobs/batches/batch_1.json`
   - 系统会解析并列出 20 个作业，显示每个作业的蛋白序列与修饰
5. **核对**（重要）：上传后浏览器会显示每个 job 的序列和 `modifications`。
   抽查 1–2 个，确认：
   - `cdk2_phos` 的 modification 是 `CCD_TPO` @ `160`
   - `mapk1_phos` 是 `CCD_TPO@185` + `CCD_PTR@187`
   - `cdk7_phos` 是 `CCD_TPO` @ `170`
6. **提交**：点击「Submit / Confirm」
7. **等待**：状态从 `Running` → `Complete`，每个作业通常几分钟到半小时
8. **下载**：完成后页面提供批量下载，得到 `folds_YYYY_MM_DD_HH_MM.zip`

## 三、下载文件放哪里

把下载的 zip **原样**放到：

```
results/af3/downloads/
```

文件名保持 `folds_*.zip` 格式（脚本按此前缀识别）。

## 四、回来后我这边做什么

你放好 zip 后，我执行：

```bash
# 1. 解析新结果（幂等，可重复跑）
python3 scripts/af3/download_af3_results.py --extract

# 2. 重跑位移比较（成对从 3 组 → 6 组）
python3 scripts/analysis/compare_af3_vs_pdb.py --all-models
python3 scripts/analysis/plot_af3_vs_pdb.py

# 3. 重跑置信度分析
python3 scripts/analysis/analyze_af3_ptm_effect.py
```

## 五、常见问题

**Q: 提示超配额怎么办？**
等第二天再传下一批。配额按 UTC 日期重置。

**Q: 上传报格式错误？**
检查 JSON 是否为数组（`[` 开头）、`modelSeeds` 是否为字符串数组（如 `["1"]`）。
当前文件已与历史成功提交的格式逐字段比对一致，正常不会报错。

**Q: 一个 job 生成几个模型？**
5 个（model_0 ~ model_4）。我们的分析脚本会遍历全部模型取平均。

**Q: 能不能改成一次传 60 个？**
可以试，但若被配额截断，需自行判断哪些已提交。建议按批次传，便于追踪。

**Q: 作业名大小写？**
服务端会把作业名转小写（`CDK2_phos_seed1` → `cdk2_phos_seed1`）。
解析脚本已做小写归一化，不影响匹配。

## 六、本批作业对照表

| 蛋白 | 基因 | 序列长 | 修饰位点(UniProt) | 类型 |
|---|---|---|---|---|
| CDK2 | cdk2 | 298 | T160 | 成对 |
| ERK2 | mapk1 | 360 | T185 + Y187 | 成对 |
| Aurora-A | aurka | 403 | T288 | 成对 |
| IRF3 | irf3 | 427 | S386 + S402 | 成对 |
| p38α | mapk14 | 360 | T180 + Y182 | 成对 |
| CDK7 | cdk7 | 346 | T170 | 成对 |
| c-Src | src | 536 | Y419 | 单结构 |
| AMPKα2 | prkaa2 | 552 | T172 | 单结构 |
| PKA-Cα | prkaca | 351 | T198 | 单结构 |
| CTD | polr2a | 1970 | S1614 + S1621 | 单结构 |

每个蛋白 6 个作业：phos × 3 种子 + unphos × 3 种子。
