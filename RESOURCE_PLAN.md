# PTM Structure Effect Research — 资源分配方案

## 总体架构

```
┌─────────────────────────────────────────────────────────┐
│                   MacBook (本地)                         │
│  • 数据采集 & 预处理                                     │
│  • AF3 作业提交 (Web API)                                │
│  • 结果下载 & 评估指标计算                                │
│  • 统计分析 & 可视化                                     │
│  • 轻量 MD 轨迹分析 (< 10GB)                             │
│  • 代码开发 & 报告撰写                                   │
└──────────────────────┬──────────────────────────────────┘
                       │ SSH/SCP
┌──────────────────────┴──────────────────────────────────┐
│              AutoDL (云端 GPU)                           │
│  • MD 模拟生产运行 (GROMACS)                             │
│  • ESMFold 推理                                          │
│  • 增强采样 (GaMD/Metadynamics)                          │
│  • 大规模轨迹分析 & PCA                                  │
│  • Markov 状态模型构建                                   │
└─────────────────────────────────────────────────────────┘
```

---

## AutoDL 云端配置需求

### 推荐租用方案

| 用途 | 推荐 GPU | VRAM | CPU | RAM | 存储 | 预计时长 | AutoDL 参考价 |
|------|---------|------|-----|-----|------|---------|--------------|
| **MD生产运行** | RTX 4090 × 1 | 24 GB | 8 vCPU | 32 GB | 200 GB | 持续租用 | ~¥2.5/h |
| **ESMFold推理** | A100 40G × 1 | 40 GB | 12 vCPU | 64 GB | 100 GB | 按需 (~2h) | ~¥6/h |
| **增强采样** | RTX 4090 × 2 | 24 GB×2 | 16 vCPU | 64 GB | 200 GB | 按需 | ~¥5/h |
| **大轨迹分析** | RTX 4090 × 1 | 24 GB | 16 vCPU | 64 GB | 500 GB | 按需 | ~¥2.5/h |

### 镜像选择

推荐 AutoDL 官方镜像：**`PyTorch 2.3.0 + CUDA 12.1 + Python 3.10 + GROMACS 2024`**

自安装补充：
```bash
conda install -c conda-forge gromacs=2024 mdtraj mdanalysis
conda install -c bioconda biopython
pip install fair-esm  # ESMFold
pip install pyemma deeptime  # MSM 分析
```

### 预算估算

| 任务 | 实例类型 | 时长 (h) | 时价 (¥) | 小计 (¥) |
|------|---------|---------|---------|---------|
| MD 批量生产 (4×500ns) | 4090 | ~200 | 2.5 | 500 |
| MD 增强采样 (2×500ns) | 4090×2 | ~100 | 5 | 500 |
| ESMFold 推理 | A100 | ~3 | 6 | 18 |
| 轨迹分析 | 4090 | ~40 | 2.5 | 100 |
| 缓冲余量 | — | — | — | 200 |
| **合计** | | | | **~¥1,318** |

> ⚠️ 可优化：多卡并行可将MD时间缩短至1/3 (~¥400-600 总费用)

---

## 详细任务分配

### 🖥️ MacBook 本地任务

| 任务 | 工具/库 | 预计时长 |
|------|--------|---------|
| 1. PDB 批量下载 | `requests` + BioPython | 30 min |
| 2. FASTA 序列提取 | BioPython `PDB.parse` | 10 min |
| 3. AF3 JSON 批量生成 | Python 脚本 | 1 h |
| 4. AF3 结果下载 & 解析 | Python + JSON | 2 h |
| 5. RMSD/TM-score 计算 | MDAnalysis + TMscore | 2 h |
| 6. 统计分析 & 可视化 | pandas + matplotlib/seaborn | 4 h |
| 7. 磷酸化分类器训练 | scikit-learn + XGBoost | 2 h |
| 8. 报告撰写 | — | 8 h |

### ☁️ AutoDL 云端任务

| 任务 | 工具 | GPU | 预计时长 |
|------|------|-----|---------|
| 1. ESMFold 预测 (45蛋白) | ESMFold v1 | A100 40G | 2-3 h |
| 2. MD 非磷酸化生产 (4体系×500ns) | GROMACS | 4090 | 60-80 h |
| 3. MD 磷酸化生产 (4体系×500ns) | GROMACS | 4090 | 60-80 h |
| 4. MD 磷酸酶模拟 (4体系×500ns) | GROMACS | 4090 | 60-80 h |
| 5. 增强采样 GaMD (2体系) | GROMACS+GaMD | 4090×2 | 50-70 h |
| 6. 轨迹精简 (stride=10) | GROMACS trjconv | 4090 | 12 h |
| 7. PCA & FEL 计算 | GROMACS + Python | 4090 | 8 h |
| 8. MSM 构建 | PyEMMA | 4090 | 6 h |

---

## 执行策略：三阶段串并行

```
阶段A (Week 1-2): 本地为主
  ├── 数据集构建 ═══════════ MacBook
  ├── AF3 JSON 批量脚本 ════ MacBook
  ├── AF3 作业提交 ═════════ MacBook → AF3 Server
  └── [等待AF3结果期间] ESMFold ═ AutoDL A100

阶段B (Week 2-5): 云端为主
  ├── MD 模拟 (phos) ═══════ AutoDL 4090 #1
  ├── MD 模拟 (unphos) ═════ AutoDL 4090 #2
  ├── 增强采样 ═════════════ AutoDL 4090×2
  └── AF3 结果分析 ═══════════ MacBook (并行)

阶段C (Week 5-8): 混合
  ├── 轨迹分析 ═════════════ AutoDL 4090 + MacBook
  ├── 综合评估 ═════════════ MacBook
  └── 报告撰写 ═════════════ MacBook
```
