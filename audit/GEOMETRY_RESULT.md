# 几何分解结果：AF3 的位移预测与真实变化零相关

日期：2026-10-27　脚本：`scripts/analysis/geometry_decomposition.py`

## 一、方法

v1 尝试用"方向一致性 |Σd|/Σ|d|"衡量位移类型，但全残基 Kabsch 叠合会令质心对齐、
Σd ≡ 0，指标恒为零（数学上的必然）。**v2 改用分块叠合**：

- `r_full`：全残基一次叠合后的 RMSD（整体形变）
- `r_block`：序列切 3 块、各自叠合后的块内 RMSD 均值（局部柔性）
- `hinge = 1 − r_block/r_full`：块间相对运动占比
  - hinge 高 → 差异主要来自块与块的相对位置（铰链/域运动）
  - hinge 低 → 差异主要是块内形变（局部重排）

## 二、核心发现：预测与实验零相关

| 检验 | 结果 |
|---|---|
| Pearson r（预测位移 vs 实验位移） | **+0.0016** |
| 线性斜率 | **0.0007**（实验每增 1 Å，AF3 仅增 0.001 Å） |
| Spearman 秩相关 | +0.064（p = 0.820） |
| 剔除 IRF3 后 | r = +0.134（p = 0.649），仍不相关 |

**AF3 的预测位移近似常数**：均值 **0.49 ± 0.33 Å**，而实验位移跨度 0.25–2.37 Å。

| 蛋白 | 实验位移 | AF3 预测 |
|---|---|---|
| Aurora-A | 2.37 Å | 0.28 Å |
| c-Src | 2.36 Å | 0.25 Å |
| AMPKα2 | 2.22 Å | 0.17 Å |
| … | … | … |
| PAK1 | 0.25 Å | 0.19 Å |

实验位移相差 9 倍的三个蛋白，AF3 给出的预测值几乎相同（0.17–0.28 Å）。

## 三、这修正了此前的表述

| 旧表述 | 新表述 |
|---|---|
| "AF3 系统性低估幅度 69%" | **"AF3 的位移幅度预测与真实变化零相关"** |
| 隐含"AF3 知道该动多少但保守" | 实际是"AF3 无论真实变化多大，都给约 0.5 Å" |

"低估 69%"是 |pred/exp| 在 exp 变化时的数学必然结果，不是 AF3 的行为特征。
本质问题是**预测值不随真实变化的幅度而变**。

## 四、副发现：铰链运动比局部形变更难预测

按 hinge 指数二分：

| 组 | n | 幅度比（中位） | 幅度比（均值） |
|---|---|---|---|
| hinge 高（铰链主导） | 8 | **0.29** | 0.31 |
| hinge 低（局部形变主导） | 7 | **0.62** | 0.83 |

**相差约 2.7 倍**，方向符合"AF3 擅长局部重排、不擅长域间运动"的预期。

但统计不显著：Welch t p = 0.171；Mann-Whitney p = 0.189。
效应量大（Cohen's d ≈ 0.88）但样本量不足，且 hinge 低组内含 IRF3（幅度比 2.70）这个离群点。

连续相关：r = −0.434（p = 0.106）；剔除 IRF3 后 r = −0.421（p = 0.134），稳定但不显著。

## 五、结论的可靠性分级

| 结论 | 强度 | 依据 |
|---|---|---|
| AF3 位移预测与真实变化零相关 | **强** | r=0.002，n=15，剔除离群点仍不显著，无需统计检验即可从散点图看出 |
| 铰链运动更难预测 | **中（趋势）** | 2.7 倍差异 + 稳定负相关，但 p>0.1，需扩样验证 |

## 六、对论文写作的要求

✅ 可以写：
- "AF3-predicted displacement magnitudes were uncorrelated with the experimentally
   observed shifts (r = 0.002, n = 15), remaining near 0.5 Å regardless of the
   true magnitude."
- "We observed a trend toward poorer performance on hinge-type motions
   (amplitude ratio 0.29 vs 0.62), though this did not reach significance
   (p = 0.17)."

❌ 不可写：
- "hinge-type motions are significantly harder for AF3"
- 任何以"低估 X%"为主语的表述

## 七、产物

| 文件 | 内容 |
|---|---|
| `results/figures/Fig6_core_result.pdf` | 双面板：零相关图 + hinge 效应 |
| `results/analysis/geometry_metrics.json` | 逐蛋白几何指标 |
| `scripts/analysis/geometry_decomposition.py` | 可复现脚本（含 v1→v2 迭代日志） |
