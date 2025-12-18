#!/usr/bin/env python3
"""
Step 8: Phyvalue 统计关联分析 (Mestre et al., 2020 方法)

基于RT-蛋白簇共现矩阵进行统计分析，筛选与Retron显著关联的蛋白簇。

核心算法 (基于Mestre et al., 2020 Supplementary Figure S3):
1. 计算每个蛋白簇在RT邻近出现的频率
2. 通过随机置换检验 (Permutation test) 计算背景分布
3. 计算Phyvalue = 观察值 / 背景期望值
4. 筛选Phyvalue > 阈值的显著关联蛋白簇

简化版方法 (无进化树时):
- 使用共现频率作为主要指标
- 通过出现在多个不同基因组中来验证非随机性

输入: Step7的rt_cluster_matrix.tsv和cluster_stats.tsv
输出:
  - significant_clusters.tsv: 显著关联的蛋白簇
  - association_scores.tsv: 所有蛋白簇的关联评分
  - phyvalue_analysis.tsv: Phyvalue分析结果
"""

import argparse
import pandas as pd
import numpy as np
import sys
import logging
from pathlib import Path
from collections import Counter, defaultdict
import random


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


def load_matrix(matrix_file, logger):
    """加载共现矩阵"""
    matrix_df = pd.read_csv(matrix_file, sep='\t', index_col=0)
    logger.info(f"加载共现矩阵: {matrix_df.shape[0]} RTs × {matrix_df.shape[1]} 蛋白簇")
    return matrix_df


def calculate_cooccurrence_frequency(matrix_df, logger):
    """
    计算每个蛋白簇的共现频率

    返回:
        freq_df: 包含每个簇的出现统计
    """
    total_rts = len(matrix_df)
    stats = []

    for cluster in matrix_df.columns:
        occurrence = matrix_df[cluster].sum()
        freq = occurrence / total_rts

        # 计算在多少个不同基因组中出现
        rt_ids = matrix_df.index[matrix_df[cluster] == 1].tolist()
        genomes = set()
        for rt_id in rt_ids:
            # RT ID格式: genome_contig_start-end
            parts = rt_id.rsplit('_', 2)
            if len(parts) >= 1:
                genomes.add(parts[0])

        stats.append({
            'cluster': cluster,
            'occurrence_count': int(occurrence),
            'occurrence_frequency': round(freq, 4),
            'genome_count': len(genomes),
            'rt_ids': ';'.join(rt_ids[:10]) if len(rt_ids) <= 10 else f"{';'.join(rt_ids[:10])}..."
        })

    freq_df = pd.DataFrame(stats)
    freq_df = freq_df.sort_values('occurrence_count', ascending=False)

    logger.info(f"计算了 {len(freq_df)} 个蛋白簇的共现频率")

    return freq_df


def permutation_test(matrix_df, n_permutations=1000, logger=None):
    """
    随机置换检验，计算背景分布

    通过随机打乱RT-簇对应关系，计算每个簇的期望出现次数

    参数:
        matrix_df: 共现矩阵
        n_permutations: 置换次数

    返回:
        background_stats: {cluster: {'mean': float, 'std': float}}
    """
    if logger:
        logger.info(f"进行 {n_permutations} 次随机置换检验...")

    clusters = matrix_df.columns.tolist()
    n_rts = len(matrix_df)

    # 存储每次置换的结果
    permutation_counts = {cluster: [] for cluster in clusters}

    for i in range(n_permutations):
        if logger and (i + 1) % 100 == 0:
            logger.debug(f"置换进度: {i + 1}/{n_permutations}")

        # 随机打乱每一列
        for cluster in clusters:
            shuffled = matrix_df[cluster].values.copy()
            np.random.shuffle(shuffled)
            permutation_counts[cluster].append(shuffled.sum())

    # 计算背景统计
    background_stats = {}
    for cluster in clusters:
        counts = permutation_counts[cluster]
        background_stats[cluster] = {
            'mean': np.mean(counts),
            'std': np.std(counts),
            'p95': np.percentile(counts, 95),
            'p99': np.percentile(counts, 99)
        }

    return background_stats


def calculate_phyvalue(freq_df, background_stats, matrix_df, logger):
    """
    计算Phyvalue (基于Mestre et al., 2020)

    Phyvalue = 观察出现次数 / 期望出现次数

    筛选条件 (根据论文):
    1. Phyvalue > 2 (观察值是期望的2倍以上)
    2. 或者 z-score > 4 (超出背景4个标准差)

    返回:
        phyvalue_df: 包含Phyvalue分析结果
    """
    results = []

    for _, row in freq_df.iterrows():
        cluster = row['cluster']
        observed = row['occurrence_count']

        bg = background_stats.get(cluster, {'mean': 0, 'std': 1, 'p95': 0, 'p99': 0})

        # 计算Phyvalue
        if bg['mean'] > 0:
            phyvalue = observed / bg['mean']
        else:
            phyvalue = float('inf') if observed > 0 else 0

        # 计算z-score
        if bg['std'] > 0:
            z_score = (observed - bg['mean']) / bg['std']
        else:
            z_score = float('inf') if observed > bg['mean'] else 0

        # 判断是否显著
        is_significant = (phyvalue > 2) or (z_score > 4)

        # 额外检查: 是否在多个基因组中出现 (排除物种特异性)
        multi_genome = row['genome_count'] > 1

        results.append({
            'cluster': cluster,
            'observed_count': observed,
            'expected_count': round(bg['mean'], 2),
            'background_std': round(bg['std'], 2),
            'phyvalue': round(phyvalue, 2),
            'z_score': round(z_score, 2),
            'p95_threshold': round(bg['p95'], 2),
            'p99_threshold': round(bg['p99'], 2),
            'is_significant': is_significant,
            'multi_genome': multi_genome,
            'genome_count': row['genome_count'],
            'occurrence_frequency': row['occurrence_frequency']
        })

    phyvalue_df = pd.DataFrame(results)
    phyvalue_df = phyvalue_df.sort_values('phyvalue', ascending=False)

    # 统计
    n_significant = phyvalue_df['is_significant'].sum()
    n_multi_genome_sig = phyvalue_df[phyvalue_df['is_significant'] & phyvalue_df['multi_genome']].shape[0]

    logger.info(f"\nPhyvalue分析结果:")
    logger.info(f"  显著关联簇 (Phyvalue>2 or z>4): {n_significant}")
    logger.info(f"  多基因组显著簇: {n_multi_genome_sig}")

    return phyvalue_df


def filter_significant_clusters(phyvalue_df, min_occurrence=5, min_phyvalue=2.0,
                                 require_multi_genome=True, logger=None):
    """
    筛选显著关联的蛋白簇

    筛选条件:
    1. 出现次数 >= min_occurrence
    2. Phyvalue >= min_phyvalue
    3. (可选) 在多个基因组中出现
    """
    filtered = phyvalue_df.copy()

    # 应用筛选条件
    filtered = filtered[filtered['observed_count'] >= min_occurrence]
    filtered = filtered[filtered['phyvalue'] >= min_phyvalue]

    if require_multi_genome:
        filtered = filtered[filtered['multi_genome'] == True]

    # 按Phyvalue排序
    filtered = filtered.sort_values('phyvalue', ascending=False)

    if logger:
        logger.info(f"\n筛选后的显著关联簇:")
        logger.info(f"  筛选条件: 出现>=   {min_occurrence}, Phyvalue>={min_phyvalue}, 多基因组={require_multi_genome}")
        logger.info(f"  筛选结果: {len(filtered)} 个蛋白簇")

    return filtered


def calculate_association_score(phyvalue_df, logger):
    """
    计算综合关联评分

    评分维度:
    - Phyvalue权重: 40%
    - 出现频率权重: 30%
    - 基因组多样性权重: 30%
    """
    df = phyvalue_df.copy()

    # 归一化各维度
    max_phyvalue = df['phyvalue'].replace([np.inf, -np.inf], np.nan).max()
    max_freq = df['occurrence_frequency'].max()
    max_genome = df['genome_count'].max()

    if pd.isna(max_phyvalue) or max_phyvalue == 0:
        max_phyvalue = 1
    if max_freq == 0:
        max_freq = 1
    if max_genome == 0:
        max_genome = 1

    # 处理无穷大
    df['phyvalue_norm'] = df['phyvalue'].apply(
        lambda x: 1.0 if x == float('inf') else min(x / max_phyvalue, 1.0)
    )
    df['freq_norm'] = df['occurrence_frequency'] / max_freq
    df['genome_norm'] = df['genome_count'] / max_genome

    # 计算综合评分
    df['association_score'] = (
        df['phyvalue_norm'] * 0.4 +
        df['freq_norm'] * 0.3 +
        df['genome_norm'] * 0.3
    )

    # 转换为0-100分
    df['association_score'] = (df['association_score'] * 100).round(1)

    # 分级
    def classify_association(score):
        if score >= 70:
            return 'HIGH'
        elif score >= 40:
            return 'MEDIUM'
        elif score >= 20:
            return 'LOW'
        else:
            return 'VERY_LOW'

    df['association_level'] = df['association_score'].apply(classify_association)

    logger.info(f"\n关联评分分布:")
    level_counts = df['association_level'].value_counts()
    for level in ['HIGH', 'MEDIUM', 'LOW', 'VERY_LOW']:
        count = level_counts.get(level, 0)
        logger.info(f"  {level}: {count}")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Step 8: Phyvalue统计关联分析 (Mestre方法)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法
  python step08_phyvalue_analysis.py -i 07_clustering -o 08_association

  # 自定义筛选条件
  python step08_phyvalue_analysis.py -i 07_clustering -o 08_association \\
    --min-occurrence 10 --min-phyvalue 2.5

  # 完整置换检验 (更准确但更慢)
  python step08_phyvalue_analysis.py -i 07_clustering -o 08_association \\
    --n-permutations 10000

方法说明:
  本步骤采用Mestre et al., 2020的Phyvalue方法:
  1. 计算每个蛋白簇在RT邻近的出现频率
  2. 通过随机置换检验估计背景分布
  3. 计算Phyvalue = 观察值 / 期望值
  4. 筛选Phyvalue > 2的显著关联蛋白簇

  简化版: 如果没有RT进化树，使用全局置换代替基于树的移动平均
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Step7输出目录 (包含rt_cluster_matrix.tsv)')
    parser.add_argument('-o', '--output', default='08_association',
                        help='输出目录 (默认: 08_association)')

    # 统计参数
    parser.add_argument('--n-permutations', type=int, default=1000,
                        help='置换检验次数 (默认: 1000)')
    parser.add_argument('--min-occurrence', type=int, default=5,
                        help='最小出现次数阈值 (默认: 5)')
    parser.add_argument('--min-phyvalue', type=float, default=2.0,
                        help='最小Phyvalue阈值 (默认: 2.0)')
    parser.add_argument('--no-multi-genome-filter', action='store_true',
                        help='不要求在多个基因组中出现')

    parser.add_argument('--seed', type=int, default=42,
                        help='随机数种子 (默认: 42)')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    # 设置随机种子
    np.random.seed(args.seed)
    random.seed(args.seed)

    # 验证输入
    input_dir = Path(args.input)
    matrix_file = input_dir / "rt_cluster_matrix.tsv"

    if not matrix_file.exists():
        logger.error(f"共现矩阵文件不存在: {matrix_file}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 8: Phyvalue统计关联分析 (Mestre方法)")
    logger.info("=" * 60)
    logger.info(f"输入目录: {input_dir}")
    logger.info(f"输出目录: {output_dir}")
    logger.info(f"置换次数: {args.n_permutations}")
    logger.info(f"最小出现次数: {args.min_occurrence}")
    logger.info(f"最小Phyvalue: {args.min_phyvalue}")

    # 加载共现矩阵
    matrix_df = load_matrix(matrix_file, logger)

    # 计算共现频率
    logger.info("\n计算共现频率...")
    freq_df = calculate_cooccurrence_frequency(matrix_df, logger)

    freq_file = output_dir / "cooccurrence_frequency.tsv"
    freq_df.to_csv(freq_file, sep='\t', index=False)
    logger.info(f"✓ 共现频率: {freq_file}")

    # 随机置换检验
    logger.info("\n进行随机置换检验...")
    background_stats = permutation_test(matrix_df, args.n_permutations, logger)

    # 计算Phyvalue
    logger.info("\n计算Phyvalue...")
    phyvalue_df = calculate_phyvalue(freq_df, background_stats, matrix_df, logger)

    phyvalue_file = output_dir / "phyvalue_analysis.tsv"
    phyvalue_df.to_csv(phyvalue_file, sep='\t', index=False)
    logger.info(f"✓ Phyvalue分析: {phyvalue_file}")

    # 筛选显著关联簇
    logger.info("\n筛选显著关联蛋白簇...")
    significant_df = filter_significant_clusters(
        phyvalue_df,
        min_occurrence=args.min_occurrence,
        min_phyvalue=args.min_phyvalue,
        require_multi_genome=not args.no_multi_genome_filter,
        logger=logger
    )

    if not significant_df.empty:
        sig_file = output_dir / "significant_clusters.tsv"
        significant_df.to_csv(sig_file, sep='\t', index=False)
        logger.info(f"✓ 显著关联簇: {sig_file} ({len(significant_df)} 个)")

    # 计算综合关联评分
    logger.info("\n计算综合关联评分...")
    scored_df = calculate_association_score(phyvalue_df, logger)

    scores_file = output_dir / "association_scores.tsv"
    scored_df.to_csv(scores_file, sep='\t', index=False)
    logger.info(f"✓ 关联评分: {scores_file}")

    # 输出高关联蛋白簇
    high_assoc = scored_df[scored_df['association_level'].isin(['HIGH', 'MEDIUM'])]
    if not high_assoc.empty:
        high_file = output_dir / "high_association_clusters.tsv"
        high_assoc.to_csv(high_file, sep='\t', index=False)
        logger.info(f"✓ 高关联簇: {high_file} ({len(high_assoc)} 个)")

    # 统计报告
    stats_file = output_dir / "analysis_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("Phyvalue统计关联分析报告\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"输入统计:\n")
        f.write(f"  RT数量: {matrix_df.shape[0]}\n")
        f.write(f"  蛋白簇数量: {matrix_df.shape[1]}\n\n")
        f.write(f"参数设置:\n")
        f.write(f"  置换次数: {args.n_permutations}\n")
        f.write(f"  最小出现次数: {args.min_occurrence}\n")
        f.write(f"  最小Phyvalue: {args.min_phyvalue}\n")
        f.write(f"  要求多基因组: {not args.no_multi_genome_filter}\n\n")
        f.write(f"结果统计:\n")
        f.write(f"  显著关联簇 (筛选后): {len(significant_df)}\n")
        f.write(f"  高关联簇 (HIGH/MEDIUM): {len(high_assoc)}\n\n")
        f.write(f"关联等级分布:\n")
        level_counts = scored_df['association_level'].value_counts()
        for level in ['HIGH', 'MEDIUM', 'LOW', 'VERY_LOW']:
            count = level_counts.get(level, 0)
            f.write(f"  {level}: {count}\n")

        if not significant_df.empty:
            f.write(f"\n显著关联簇Top 10:\n")
            for _, row in significant_df.head(10).iterrows():
                f.write(f"  {row['cluster']}: Phyvalue={row['phyvalue']}, "
                       f"出现={row['observed_count']}, 基因组={row['genome_count']}\n")

    logger.info(f"✓ 分析报告: {stats_file}")

    # 总结
    logger.info("\n" + "=" * 60)
    logger.info("Phyvalue统计关联分析完成!")
    logger.info(f"  分析蛋白簇数: {len(phyvalue_df)}")
    logger.info(f"  显著关联簇: {len(significant_df)}")
    logger.info(f"  高关联簇: {len(high_assoc)}")
    logger.info(f"  输出目录: {output_dir}")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
