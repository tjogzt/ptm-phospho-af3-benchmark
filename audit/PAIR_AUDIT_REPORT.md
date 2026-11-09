# 基准配对质量筛查报告（最终）

筛查脚本：`scripts/audit_pairs.py`　执行日期：2026-10-22
方法：直接解析 PDB 头部 TITLE/COMPND + MODRES 记录 + 残基名 + 链序列，四项硬核验。

## 核心结论

**18 个基准蛋白中，仅 6 个的"磷酸化结构"真实含 SEP/TPO/PTR 修饰。**

其中 9 组成对结构中只有 5 组 phos 结构为真；叠加"AF3 数据已返回"与"存在实验位移"两个条件后，
**可用于 AF3-vs-实验验证的配对仅剩 3 组**（CDK2 / ERK2 / Aurora-A）。

## 一、成对结构组（9 组）

| 蛋白 | unphos | phos | phos 真实磷酸化 | phos 实际内容 | 判定 |
|---|---|---|---|---|---|
| CDK2 | 1HCL | 1JST | ✅ 1×TPO | PHOSPHORYLATED CDK2-CYCLIN A | **PASS** |
| ERK2 | 1ERK | 2ERK | ✅ 2× | PHOSPHORYLATED MAP KINASE ERK2 | **PASS** |
| Aurora-A | 1MQ4 | 1OL5 | ✅ 1× | AURORA-A PHOSPHORYLATED ON T287/T288 | **PASS** |
| PKA-Cα | 1CTP | 1ATP | ✅ 2× | PKA 催化亚基 | ⚠ 链 ID 为 **E/I**，非声明的 A；AF3 未返回 |
| PYGM | 1GPA | 1GPA | ✅ 1×SEP14 | 同一结构 | ⚠ 两态同源，无实验位移 |
| CheY | 3CHY | 1FQW | ❌ 0 | ACTIVATED CHEY（BeF3 模拟） | ✗ AF3 不支持 Asp-P |
| CDK7 | 1UA2 | 7ML6 | ❌ 0 | **细菌蛋白 CALU17**（非 CDK7） | ✗ **PDB ID 错误** |
| IRF3 | 1QWT | 1T2K | ❌ 0 | **DNA 结合域复合物**，N 端 DBD vs C 端 TAD | ✗ **结构域错配** |
| p38α | 1WFC | 3D83 | ❌ 0 | **p38-联苯酰胺抑制剂复合物** | ✗ **非磷酸化结构** |

## 二、单体结构组（7 组）

| 蛋白 | PDB | phos 真实磷酸化 | 实际内容 | 判定 |
|---|---|---|---|---|
| CTD (POLR2A) | 4JXT | ✅ 1× | RPRD1A CID + 磷酸化 RPB1-CTD 肽 | **PASS**（链为肽链） |
| EGFR | 3W2S | ❌ 0 | EGFR KINASE DOMAIN WITH **COMPOUND4** | ✗ 化合物复合物 |
| c-Src | 1Y57 | ❌ 0 | **"UNPHOSPHORYLATED C-SRC"**（标题明写非磷酸化） | ✗ |
| AMPKα2 | 2H6D | ❌ 0 | AMPK α2 激酶结构域 | ✗ 无磷酸化 |
| CFTR | 5UAK | ❌ 0 | **"DEPHOSPHORYLATED, ATP-FREE"**（标题明写去磷酸化） | ✗ |
| Hsp27 | 4MHR | ❌ 0 | **细菌 S. alaskensis EctD 的 apo 形式**（非 Hsp27） | ✗ **蛋白身份错误** |
| PhKγ | 1PHK | ❌ 0 | 磷酸化酶激酶催化域 + 核苷酸 | ✗ 无磷酸化 |

## 三、负对照（2 组）

| 蛋白 | PDB | 判定 |
|---|---|---|
| Lysozyme | 1LYS | PASS（人工引入位点，符合阴性对照设计） |
| Ubiquitin | 5OXL / 1UBQ | ⚠ 5OXL 为 PINK1-泛素复合物，非 76 aa 成熟泛素 |

## 四、筛选后的可用配对（3 组，15 模型平均）

| 蛋白 | 实验位移 | AF3 预测位移 | 方向余弦 | 幅度 r |
|---|---|---|---|---|
| CDK2 | 2.70 Å | 0.63 Å | +0.440 | +0.382 |
| ERK2 | 1.61 Å | 1.00 Å | +0.661 | +0.580 |
| Aurora-A | 0.80 Å | 0.27 Å | +0.471 | +0.805 |

**初步结论**：AF3 能捕捉磷酸化诱导位移的**方向**（余弦 0.44–0.66），但 3 组中有 2 组
显著**低估幅度**。样本量过小（n=3），不足以支撑统计推断。

## 五、处置建议（按优先级）

1. **立即剔除**：CDK7(7ML6)、IRF3(1T2K)、Hsp27(4MHR)、p38α(3D83)、c-Src(1Y57)、CFTR(5UAK)、
   EGFR(3W2S)、AMPKα2(2H6D)、PhKγ(1PHK) —— 这些"磷酸化结构"不成立，基于它们的任何
   phos/unphos 差异都不能归因于磷酸化。
2. **修正链 ID**：PKA-Cα 的 1ATP/1CTP 应使用链 E（非 A）。
3. **替换 PDB**：为 CDK7、IRF3、p38α、Hsp27、c-Src、CFTR、EGFR 重新检索真实的磷酸化结构
   （建议以 PDBe 的 `pdbx_struct_mod_residue` / 磷酸化注释为检索条件，而非按蛋白名检索）。
4. **重建筛选门禁**：将 `scripts/audit_pairs.py` 作为基准入库的前置检查，任何新增配对
   必须 C1–C4 全部通过。

## 六、为什么原基准会出这些问题

原流程按"蛋白名 + 磷酸化关键词"检索 PDB 并人工填入 ID，未对候选结构做机械核验，
导致三类错误混入：
- **PDB ID 张冠李戴**（7ML6、4MHR）—— ID 指向完全不同的蛋白
- **非磷酸化结构被当作磷酸化态**（3D83、3W2S、2H6D、1Y57、5UAK、1PHK）
- **结构域错配**（IRF3 的 1QWT/1T2K）

这三类错误都无法通过"看 pLDDT 数值"发现，只能靠结构性核验。本脚本即为该门禁。
