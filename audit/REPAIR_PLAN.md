# 基准修复方案（基于真实磷酸化结构检索）

检索方式：RCSB Search API 按 UniProt accession 列出全部 entry → HTTP Range 取 PDB 头部 →
解析 MODRES 记录 → 仅保留真含 SEP/TPO/PTR 者。机器可读结果见 `audit/phos_structure_candidates.json`。

## 一、可直接替换的结构（6 个蛋白）

| 蛋白 | 原 phos PDB（错误） | 替代候选 | 分辨率 | 关键磷酸化位点 | 证据 |
|---|---|---|---|---|---|
| **IRF3** | 1T2K（DNA 结合域，非磷酸化） | **3A77** | 1.8 Å | SEP@386, 402 | 标题："THE CRYSTAL STRUCTURE OF **PHOSPHORYLATED** IRF-3" |
| | | 7JFL | 1.68 Å | SEP@386 | "PHOSPHORYLATED IRF-3 BOUND TO CBP" |
| **P38α** | 3D83（抑制剂复合物，非磷酸化） | **3NNX** | 2.28 Å | SEP@182 | "CRYSTAL STRUCTURE OF **PHOSPHORYLATED** P38 ALPHA" |
| | | 4EWQ | 2.1 Å | SEP@180 | p38α + 哒嗪类抑制剂（含 pT180） |
| **c-Src** | 1Y57（标题明写 UNPHOSPHORYLATED） | **1YI6** | 2.0 Å | SEP@335, 416 | "C-TERM TAIL SEGMENT OF HUMAN TYROSINE KINASE (258-533)" |
| **AMPKα2** | 2H6D（无磷酸化） | **6B1U** | 2.77 Å | SEP@172 | "FULL-LENGTH HUMAN AMPK (A2B1G1)" |
| **CDK7** | 7ML6（细菌蛋白 CALU17） | **8PYR** | 2.15 Å | SEP@170 | 含 SEP@164/170 |
| **PKA-Cα** | 1ATP（链 ID 应为 E） | **6QJ7** | 1.69 Å | SEP@197 | PKA RIA 全酶；1ATP 保留但改用链 E |

## 二、无法修复的（2 个蛋白）

| 蛋白 | 检索结果 | 建议处置 |
|---|---|---|
| **Hsp27** (HSPB1) | 7 个 entry 中 **0 个**含磷酸化修饰 | 移出磷酸化结构基准；或降级为"无序蛋白阴性对照"（其 pLDDT 58–60 本就提示高柔性） |
| **PhKγ** (PHKG1) | UniProt 仅 1 个 entry，**无磷酸化结构** | 移出基准 |

## 三、需要重新定义靶标的（1 个蛋白）

**CFTR**：57 个 entry 中 8 个含磷酸化，但**全部是 R 域肽段与 14-3-3 蛋白的复合物**
（5D3E / 7QI1 / 6HEP / 5D2D / 5D3F，均为 pS753/pS768 等肽段）。
不存在全长 CFTR 的磷酸化实验结构。

→ 建议：CFTR 的 disorder→order 命题要么改为"R 域肽段 + 14-3-3 结合"体系（与 AF3 单体预测不可比），
要么放弃该条目。**当前全长 1480 aa 的 AF3 预测缺乏可对照的实验结构。**

## 四、检索中新发现的 benchmark 错误（补遗）

**CDK7 的 unphos 参照 1UA2 本身含 SEP@170**——即原 benchmark 用作"非磷酸化"参照的结构
实际上已是磷酸化状态。因此 CDK7 配对（1UA2 vs 7ML6）**两端都错**：
- unphos 端已是磷酸化（1UA2 含 pT170）
- phos 端是完全无关的细菌蛋白（7ML6）

## 五、修复后的基准构成（预期）

| 类别 | 修复前 | 修复后 |
|---|---|---|
| 成对结构（两态均真实） | 9 组（实际可用 3 组） | **6 组**：CDK2, ERK2, Aurora-A, PKA(改链E), IRF3(→3A77), P38α(→3NNX) |
| 单结构磷酸化 | 7 组（仅 1 组真实） | **3 组**：c-Src(→1YI6), AMPKα2(→6B1U), CDK7(→8PYR) |
| 待移出 | — | Hsp27, PhKγ, CFTR（无合适实验结构） |
| 阴性对照 | 2 组 | 2 组（Lysozyme, Ubiquitin；UBB 需改 76 aa） |

## 六、执行顺序建议

1. **下载替代结构**并跑 `scripts/audit_pairs.py` 复核（确认 C1–C4 全过）
2. **重建 unphos 参照**：CDK7 需另找非磷酸化结构（1UA2 已含 pT170）
3. **重新生成 AF3 作业**（`generate_af3_jobs.py` 需改为读 `data/curated/` 的 UniProt 序列）
4. **提交并重跑**，然后重做 Fig2 的位移比较
