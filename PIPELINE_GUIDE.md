# Retron RT 挖掘流程 v2.1 完整指南

## 流程概述

本流程用于从细菌基因组中挖掘 Retron 系统，基于 **Mestre et al., 2020 NAR** 方法学，采用统计关联分析 (Phyvalue) 识别 Retron 相关蛋白。

### 流程架构图

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        Retron RT 挖掘流程 v2.1                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌──────────────────┐      ┌──────────────────┐                            │
│  │ 方案A: DIAMOND   │  OR  │ 方案B: HMMSEARCH │  ← Step 1: RT 搜索         │
│  │ (序列相似性搜索)  │      │ (HMM模型搜索)    │                            │
│  └────────┬─────────┘      └────────┬─────────┘                            │
│           │                         │                                       │
│           └───────────┬─────────────┘                                       │
│                       ▼                                                     │
│           ┌──────────────────────┐                                         │
│           │ Step 2: Tier 分层分类 │  ← NAXXH/VTG/YADD motif 检测            │
│           │ Tier1: 高置信度       │                                         │
│           │ Tier2: 潜在新型       │                                         │
│           │ Tier3: 排除候选       │                                         │
│           └──────────┬───────────┘                                         │
│                      ▼                                                      │
│           ┌──────────────────────┐                                         │
│           │ Step 3: 系统发育分析  │  ← MAFFT L-INS-i + trimal + FastTree    │
│           │ 距离标注置信度        │                                         │
│           └──────────┬───────────┘                                         │
│                      ▼                                                      │
│           ┌──────────────────────┐                                         │
│           │ Step 4: 邻近蛋白提取  │  ← Mestre 方法: ±20-30kb               │
│           └──────────┬───────────┘                                         │
│                      ▼                                                      │
│           ┌──────────────────────┐                                         │
│           │ Step 5: MMseqs2 聚类  │  ← 30% identity, 80% coverage          │
│           └──────────┬───────────┘                                         │
│                      ▼                                                      │
│           ┌──────────────────────┐                                         │
│           │ Step 6: Phyvalue 分析 │  ← 统计关联分析                         │
│           └──────────┬───────────┘                                         │
│                      ▼                                                      │
│           ┌──────────────────────┐                                         │
│           │ Step 7: 类型分类报告  │  ← Retron Type I-XIII                   │
│           └──────────┬───────────┘                                         │
│                      ▼                                                      │
│           ┌──────────────────────┐                                         │
│           │ Step 8: ncRNA 验证    │  ← 可选: 共变分析                       │
│           └──────────────────────┘                                         │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Step 0: 构建 Retron RT HMM 模型 (可选但推荐)

### 目的
从12条已知 Retron RT 参考序列构建特异性 HMM 模型，用于更精确的 RT 搜索。

### 优势
- 相比通用 Pfam RVT_1 模型，Retron RT 特异 HMM 能更准确识别 NAXXH/VTG motif 位置
- 提供 `env_from/env_to` 边界，精确定位 RT 核心区域

### 参考序列 (12条)
| 序列名 | 来源物种 | 说明 |
|--------|----------|------|
| Ec67 | E. coli | 经典 Retron-Ec67 |
| Ec86 | E. coli | 经典 Retron-Ec86 |
| Ec73 | E. coli | |
| Ec83 | E. coli | |
| Ec48 | E. coli | |
| Ec107 | E. coli | |
| Ec78 | E. coli | |
| MX65 | Myxococcus xanthus | 粘细菌来源 |
| ML162 | Melittangium | 粘细菌来源 |
| Sa163 | Stigmatella | |
| Vc95 | Vibrio | |
| Se72 | Salmonella | |

### 运行
```bash
cd /home/teng/claude_code/retron/pipeline_v2/retron_scripts
bash build_retron_hmm.sh
```

### 输出
```
/home/teng/claude_code/retron/database/retron_hmm/
├── reference_rt.fasta           # 合并的参考序列
├── reference_rt_aligned.fasta   # MAFFT L-INS-i 比对结果
├── reference_rt_aligned.sto     # Stockholm 格式
├── Retron_RT.hmm               # HMM 模型文件
├── Retron_RT.hmm.h3f           # hmmsearch 索引文件
├── Retron_RT.hmm.h3i
├── Retron_RT.hmm.h3m
└── Retron_RT.hmm.h3p
```

---

## Step 1: RT 搜索

### 方案选择

| 方案 | 脚本 | 方法 | 优势 | 适用场景 |
|------|------|------|------|----------|
| A | `step01_search.py` | DIAMOND/BLAST | 速度快，批处理 | 大规模初筛 |
| B | `step01_hmmsearch.py` | HMMER (自建HMM) | 精确识别核心区域 | 精确挖掘 |

### 方案 A: DIAMOND 搜索
```bash
python step01_search.py \
    -q /home/teng/claude_code/retron/database/query_fasta \
    -d /home/teng/claude_code/retron/faa_failes \
    -o 01_search \
    --min-identity 25 \
    --min-coverage 40 \
    -t 8
```

### 方案 B: HMM 搜索 (推荐)
```bash
# 先构建 HMM 模型 (如果还没构建)
bash build_retron_hmm.sh

# 运行 HMM 搜索
python step01_hmmsearch.py \
    -m /home/teng/claude_code/retron/database/retron_hmm/Retron_RT.hmm \
    -d /home/teng/claude_code/retron/faa_failes \
    -o 01_hmm_search \
    --min-score 25 \
    --evalue 1e-5 \
    -t 8
```

### 输出
```
01_search/reports/complete_results.tsv
```

关键列说明 (HMM搜索版本):
- `env_from`, `env_to`: RT 核心区域边界 (氨基酸位置)
- `hmm_score`: HMM 匹配分数
- `hmm_coverage`: HMM 模型覆盖率

---

## Step 2: Tier 分层分类

### 目的
基于保守 motif 对 RT 候选进行**软分类**（不进行硬过滤），保留潜在新型 Retron。

### Motif 检测

| Motif | 模式 | 意义 |
|-------|------|------|
| NAXXH | `[NQDE]A..[HY]` | Retron RT 特征 (Region X) |
| VTG | `V[TS]G` (C端150aa内) | Retron RT 特征 (Region Y) |
| YADD | `Y[AV]DD` | Group II Intron RT 标志 (排除标记) |

### Tier 分类标准

| Tier | 条件 | 置信度 | 建议操作 |
|------|------|--------|----------|
| **Tier 1** | NAXXH + VTG (无YADD) | 高置信度 | keep - 直接进入下一步 |
| **Tier 2** | 只有 NAXXH 或 VTG (无YADD) | 潜在新型 | review - 人工审查 |
| **Tier 3** | 有 YADD 标记 或 无任何motif | 排除 | exclude - 建议排除 |

### 运行
```bash
python step02_motif_filter.py \
    -i 01_search/reports/complete_results.tsv \
    -f /home/teng/claude_code/retron/faa_failes \
    -o 02_motif_tier
```

### 输出
```
02_motif_tier/
├── motif_tier_all.tsv       # 所有候选 (含 tier 标签)
├── tier1_candidates.tsv     # Tier 1: 高置信度
├── tier2_candidates.tsv     # Tier 2: 潜在新型
├── tier3_excluded.tsv       # Tier 3: 排除候选
├── motif_filtered.tsv       # 兼容旧流程 (默认包含所有tier)
├── high_confidence.tsv      # 兼容旧流程 (=Tier 1)
└── tier_classification_stats.txt  # 统计报告
```

### 下一步建议
- 保守策略: 只使用 `tier1_candidates.tsv`
- 探索策略: 使用 `motif_tier_all.tsv` 或 `tier1 + tier2`

---

## Step 3: 系统发育分析

### 方法改进 (v2.1)

| 组件 | 旧版 | 新版 | 改进说明 |
|------|------|------|----------|
| MAFFT | `--auto` | `--localpair --maxiterate 1000` | L-INS-i算法，适合低相似度(20-35%)序列 |
| trimal | 无 | `-automated1` | 自动修剪低质量比对区域 |
| FastTree | 默认参数 | `-lg -gamma` | LG模型 + gamma分布，更准确 |

### 置信度标注

基于到最近 Retron 参考序列的系统发育距离:

| 置信度 | 距离范围 | 建议操作 |
|--------|----------|----------|
| high | < 0.3 | keep |
| medium | 0.3 - 0.5 | keep |
| low | 0.5 - 0.75 | review |
| very_low | >= 0.75 | review |

### 运行
```bash
python step03_phylogeny.py \
    -i 02_motif_tier/tier1_candidates.tsv \
    -f /home/teng/claude_code/retron/faa_failes \
    -r /home/teng/claude_code/retron/database/rt_reference \
    -o 03_phylogeny \
    --distance-threshold 0.5 \
    --threads 8
```

### 输出
```
03_phylogeny/
├── rt_sequences_for_alignment.fasta  # 输入序列
├── rt_alignment.fasta                # MAFFT L-INS-i 比对
├── rt_alignment_trimmed.fasta        # trimal 修剪后
├── rt_tree.nwk                       # 系统发育树
├── retron_candidates.tsv             # 所有候选 (含距离标注)
├── rt_classification_summary.tsv     # 置信度汇总
└── phylogeny_stats.txt               # 统计信息
```

---

## Step 4-8: 后续分析流程

### Step 4: 邻近蛋白提取 (Mestre 方法)
```bash
python step04_extract_neighbors.py \
    -i 03_phylogeny/retron_candidates.tsv \
    -a /home/teng/claude_code/retron/antismash_files \
    -o 04_neighbors \
    --distance 20  # kb
```

### Step 5: MMseqs2 蛋白聚类
```bash
# 可能需要单独的 conda 环境
conda activate mmseqs_env

python step05_mmseqs_cluster.py \
    -i 04_neighbors \
    -o 05_clustering \
    --min-seq-id 0.3 \
    --coverage 0.8
```

### Step 6: Phyvalue 统计关联分析
```bash
python step06_phyvalue_analysis.py \
    -i 05_clustering \
    -o 06_association \
    --min-phyvalue 2.0 \
    --n-permutations 1000
```

### Step 7: 类型分类报告
```bash
python step07_final_report.py \
    -i 06_association \
    -o 07_report
```

### Step 8: ncRNA 共变验证 (可选)
```bash
python step08_ncrna_validation.py \
    -i 07_report/final_candidates.tsv \
    -f /home/teng/claude_code/retron/fasta_files \
    -o 08_ncrna \
    --phylogeny 03_phylogeny/retron_candidates.tsv
```

---

## 完整运行方案

### 方案1: 使用配置文件 (推荐)

```bash
cd /home/teng/claude_code/retron/pipeline_v2/retron_scripts

# 完整流程
python main_pipeline.py -c config.yaml

# 干运行 (预览命令)
python main_pipeline.py -c config.yaml --dry-run

# 从特定步骤开始
python main_pipeline.py -c config.yaml --start-step 3

# 跳过 ncRNA 验证
python main_pipeline.py -c config.yaml --skip-ncrna
```

### 方案2: HMM 搜索流程 (手动步骤)

```bash
cd /home/teng/claude_code/retron/pipeline_v2/retron_scripts

# Step 0: 构建 HMM 模型 (首次运行)
bash build_retron_hmm.sh

# Step 1: HMM 搜索
python step01_hmmsearch.py \
    -m /home/teng/claude_code/retron/database/retron_hmm/Retron_RT.hmm \
    -d /home/teng/claude_code/retron/faa_failes \
    -o 01_hmm_search \
    -t 8

# Step 2: Tier 分层分类
python step02_motif_filter.py \
    -i 01_hmm_search/reports/complete_results.tsv \
    -f /home/teng/claude_code/retron/faa_failes \
    -o 02_motif_tier

# Step 3: 系统发育分析
python step03_phylogeny.py \
    -i 02_motif_tier/tier1_candidates.tsv \
    -f /home/teng/claude_code/retron/faa_failes \
    -o 03_phylogeny \
    --threads 8

# Step 4-8: 继续后续分析...
```

### 方案3: 分阶段运行 (大规模数据)

```bash
# 阶段1: 搜索和过滤 (Step 1-3)
python main_pipeline.py -c config.yaml --end-step 3

# 阶段2: 聚类分析 (需要 mmseqs 环境)
conda activate mmseqs_env
python main_pipeline.py -c config.yaml --start-step 4 --end-step 5

# 阶段3: 关联分析和报告 (Step 6-8)
conda deactivate
python main_pipeline.py -c config.yaml --start-step 6
```

---

## 依赖安装

```bash
# Python 包
pip install biopython pandas pyyaml tqdm psutil

# 核心工具
conda install -c bioconda diamond mafft fasttree trimal hmmer

# MMseqs2 (可能需要单独环境)
conda create -n mmseqs_env mmseqs2
conda activate mmseqs_env

# ncRNA 验证 (可选)
conda install -c bioconda viennarna r-scape
```

---

## 输出目录结构

```
output/
├── 01_search/              # 或 01_hmm_search/
│   ├── reports/
│   │   └── complete_results.tsv
│   └── logs/
│
├── 02_motif_tier/
│   ├── motif_tier_all.tsv
│   ├── tier1_candidates.tsv
│   ├── tier2_candidates.tsv
│   ├── tier3_excluded.tsv
│   └── tier_classification_stats.txt
│
├── 03_phylogeny/
│   ├── rt_alignment.fasta
│   ├── rt_alignment_trimmed.fasta
│   ├── rt_tree.nwk
│   ├── retron_candidates.tsv
│   └── rt_classification_summary.tsv
│
├── 04_neighbors/
│   ├── all_neighbors.faa
│   └── neighbor_matrix.tsv
│
├── 05_clustering/
│   ├── rt_cluster_matrix.tsv
│   └── cluster_info.tsv
│
├── 06_association/
│   ├── phyvalue_analysis.tsv
│   └── significant_clusters.tsv
│
├── 07_report/
│   ├── final_candidates.tsv
│   └── retron_types.tsv
│
└── 08_ncrna/               # 可选
    ├── covariance_results.tsv
    └── validated_candidates.tsv
```

---

## 关键技术要点

### 1. 低相似度序列处理
- RT 酶在 20-35% 序列相似度范围内仍可能是真正的同源物
- MAFFT L-INS-i 算法专门优化低相似度比对
- 参考序列作为"坐标轴"帮助定位候选

### 2. Tier 软分类 vs 硬过滤
- **不推荐**: 直接过滤掉无 motif 的序列
- **推荐**: 使用 Tier 标签进行软分类，保留潜在新型 Retron

### 3. HMM 搜索优势
- Retron RT 特异 HMM 比通用 Pfam RVT_1 更准确
- `env_from/env_to` 提供精确的核心区域边界
- 适用于后续结构域分析

### 4. 系统发育距离解读
- 距离 < 0.3: 与已知 Retron RT 高度相似
- 距离 0.3-0.5: 中等相似，可能是同家族变体
- 距离 > 0.5: 距离较远，需要人工审查

---

## 常见问题

### Q1: DIAMOND 搜索 vs HMM 搜索 选哪个？
- **大规模初筛**: 使用 DIAMOND (速度快)
- **精确挖掘**: 使用 HMM (更准确的核心区域识别)
- **建议**: 可以两种方法并行，取并集

### Q2: Tier 2 候选如何处理？
- 人工检查序列比对
- 查看是否有部分保守的 motif 变体
- 结合 Step 3 系统发育距离判断

### Q3: trimal 未安装怎么办？
- Step 3 会自动 fallback 到未修剪的比对
- 建议安装: `conda install -c bioconda trimal`

### Q4: 如何处理大量候选 (>1000)?
- 使用批处理模式 (`--batch-size 100`)
- 分阶段运行，先完成 Step 1-3 筛选
- 只对高置信度候选进行后续分析
