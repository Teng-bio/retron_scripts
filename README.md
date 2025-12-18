# Retron系统挖掘流程 v2.0

基于Millman et al., 2020 Cell的方法论重构，整合多维度筛选策略。

## 流程概览

```
Step 1: RT搜索 (DIAMOND/BLAST)
    ↓
Step 2: Motif过滤 (NAXXH + VTG) ← 新增：关键过滤步骤
    ↓
Step 3: RT分类与排除 (Group II/DGR/CRISPR) ← 新增
    ↓
Step 4: 上游序列提取
    ↓
Step 5: RNA结构筛选 (滑动窗口 + 分支G)
    ↓
Step 6: 邻近基因提取
    ↓
Step 7: 防御岛分析 (DefenseFinder) ← 新增
    ↓
Step 8: 效应蛋白检测 ← 新增
    ↓
Step 9: 综合报告 + 置信度分级 ← 新增
```

## 快速开始

### 1. 使用配置文件运行

```bash
# 编辑配置文件
cp config_v2.yaml my_config.yaml
vim my_config.yaml

# 运行完整流程
python main_pipeline.py -c my_config.yaml

# 预览命令（不执行）
python main_pipeline.py -c my_config.yaml --dry-run
```

### 2. 命令行参数运行

```bash
# 如果已有Step 1的搜索结果
python main_pipeline.py \
  --search-results /path/to/complete_results.tsv \
  --faa-dir /path/to/faa \
  --fasta-dir /path/to/fasta \
  --antismash-dir /path/to/antismash \
  --output ./retron_v2_results

# 从头开始
python main_pipeline.py \
  --query-dir /path/to/reference_RT \
  --faa-dir /path/to/faa \
  --fasta-dir /path/to/fasta \
  --output ./retron_v2_results
```

### 3. 分步运行

```bash
# Step 2: Motif过滤
python step02_motif_filter.py \
  -i complete_results.tsv \
  -f /path/to/faa \
  -o 02_motif_filtered

# Step 3: RT分类
python step03_classify_rt.py \
  -i 02_motif_filtered/motif_filtered.tsv \
  -o 03_classified

# Step 7: 防御岛分析
python step07_defense_island.py \
  -i rna_filtered.tsv \
  -f /path/to/faa \
  -o 07_defense_island
```

## 输出文件说明

```
retron_v2_results/
├── 01_search/           # RT搜索结果
├── 02_motif_filtered/   # Motif过滤结果
│   ├── motif_filtered.tsv     # 通过过滤的RT
│   └── high_confidence.tsv    # NAXXH+VTG都有的RT
├── 03_classified/       # RT分类结果
│   ├── retron_candidates.tsv  # Retron候选
│   └── excluded_rt.tsv        # 被排除的RT
├── 04_flanking/         # 上游序列
├── 05_rna_filtered/     # RNA筛选结果
├── 06_neighbors/        # 邻近基因
├── 07_defense_island/   # 防御岛分析
├── 08_effector/         # 效应蛋白检测
├── 09_report/           # 最终报告
│   ├── final_results.tsv              # 完整结果
│   ├── high_confidence_candidates.tsv # 高置信度候选
│   └── summary_report.md              # 摘要报告
└── logs/                # 运行日志
```

## 置信度评分标准

| 维度 | 分值 | 条件 |
|------|------|------|
| Motif | 30 | NAXXH + VTG |
| Motif | 20 | NAXXH only |
| RNA茎长 | 15 | ≥20bp |
| RNA茎长 | 10 | ≥15bp |
| 分支G | 5 | 存在 |
| 效应蛋白 | 25 | 邻近有效应蛋白 |
| 防御岛 | 20 | 位于防御系统附近 |

**置信度等级:**
- HIGH (≥70分): 高置信度，可直接实验验证
- MEDIUM (45-69分): 中等置信度，建议补充分析
- LOW (25-44分): 低置信度，需更多证据
- VERY_LOW (<25分): 很可能是假阳性

## 依赖

```bash
# Python包
pip install biopython pandas pyyaml tqdm

# 外部工具
conda install -c bioconda diamond blast viennarna

# 可选：防御岛分析
conda install -c bioconda defense-finder
```

## 小样本测试建议

首次使用时建议先用小样本测试：

1. 选择50-100个基因组子集
2. 确保包含已知含Retron的菌株（如M. xanthus）
3. 运行完整流程
4. 检查高置信度候选是否包含已知Retron
5. 根据结果调整参数

## 参考文献

1. Millman A, et al. (2020) Bacterial Retrons Function In Anti-Phage Defense. Cell.
2. Zimmerly S, Wu L. (2015) An Unexplored Diversity of Reverse Transcriptases in Bacteria.
