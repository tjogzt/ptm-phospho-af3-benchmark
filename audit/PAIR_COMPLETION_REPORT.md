# 配对补齐报告（6 组 → 10 组）

日期：2026-10-24　方法：`scripts/find_unphos_structures.py` + `scripts/find_new_pairs.py`

## 核心结果

成对结构由 6 组补齐至 **10 组**，全部经五条判据核验通过。

## 一、三处"原基准误用"的反转

搜索非磷酸化结构时发现，原 benchmark 有三个结构本就是**非磷酸化态**，被错误地当成磷酸化结构使用。
把它们放回正确位置后，直接补出 3 组配对：

| 蛋白 | PDB | 原用途 | 实际内容 | 修正后用途 |
|---|---|---|---|---|
| c-Src | 1Y57 | 被当作 phos | 标题："**UNPHOSPHORYLATED** C-SRC" | **unphos 对照** ✓ |
| AMPKα2 | 2H6D | 被当作 phos | MODRES 为空，无任何磷酸化 | **unphos 对照** ✓ |
| PYGM | 1H5U | 未使用 | "GLYCOGEN PHOSPHORYLASE **B**"（非磷酸化） | **unphos 对照** ✓ |

其中 c-Src 与 AMPK 的 phos 结构另有正确替代（1YI6 / 6B1U），因此这两组配对得以成立。

## 二、新增第 10 组：CDK1

自动筛查候选激酶后选定 CDK1（与已有 CDK2、CDK7 同属 CDK 家族，可做家族内对比）：

| 状态 | PDB | 链 | 磷酸化 | 覆盖度 | 位点 |
|---|---|---|---|---|---|
| unphos | 4YC3 | A | 无 | 99% | 覆盖 T161 |
| phos | 9SKQ | A | TPO@161 | 96% | T-loop |

## 三、无法补齐的 2 个蛋白

| 蛋白 | 检索结果 | 原因 |
|---|---|---|
| PKA-Cα | 54 个 entry **全部含 pT197** | T197 磷酸化是激酶折叠与活性所必需，实验上不存在非磷酸化晶体结构 |
| CTD (POLR2A) | 58 个 entry 无合适肽段 | CTD 肽结构稀少，且多为磷酸化态或与其他蛋白复合 |

两者保留为单结构组，不参与位移验证。

## 四、10 组配对总表

| # | 蛋白 | unphos | phos | 位点(UniProt) | 效应类型 |
|---|---|---|---|---|---|
| 1 | CDK2 | 1HCL:A | 1JST:A | T160 | local_rearrangement |
| 2 | ERK2 | 1ERK:A | 2ERK:A | T185/Y187 | local_rearrangement |
| 3 | Aurora-A | 1MQ4:A | 1OL5:A | T288 | local_rearrangement |
| 4 | IRF3 | 1QWT:A | 3A77:A | S386/S402 | allosteric |
| 5 | p38α | 1WFC:A | 3PY3:A | T180/Y182 | local_rearrangement |
| 6 | CDK7 | 7B5O:J | 8PYR:A | T170 | local_rearrangement |
| 7 | c-Src | 1Y57:A | 1YI6:A | Y419 | allosteric |
| 8 | AMPKα2 | 2H6D:A | 6B1U:A | T172 | allosteric |
| 9 | PYGM | 1H5U:A | 1GPA:A | S15 | disorder_order |
| 10 | CDK1 | 4YC3:A | 9SKQ:A | T161 | local_rearrangement |

效应类型分布：local_rearrangement 6、allosteric 3、disorder_order 1。

## 五、AF3 作业状态

| 蛋白 | v2 数据 | 待补 |
|---|---|---|
| CDK2/ERK2/Aurora-A/IRF3/p38α/CDK7 | ✅ 完整（3 种 × 2 态） | — |
| c-Src / AMPKα2 | phos_seed1 已回 | 5 个 |
| PYGM / CDK1 | 无 | 各 6 个 |
| PKA / CTD（单结构） | phos_seed1 已回 | 5 个 |

待补合计 **32 个作业**，已生成 `af3_jobs/batch_4_supplement.json` 并分为 batch_4 / batch_5。

## 六、核验记录

10 组配对全部通过：phos 含 SEP/TPO/PTR、unphos 无修饰、序列覆盖 46–100%、
位点落在结构范围内、链 ID 逐一指定。详见逐条核验输出。
