#!/usr/bin/env python3
"""
Step 7: MMseqs2 聚类与共现矩阵构建 (Mestre et al., 2020 方法)

对Step6提取的所有邻近蛋白进行聚类，并构建RT-蛋白簇共现矩阵。
这是统计关联分析的基础数据结构。

方法流程:
1. 使用MMseqs2对所有邻近蛋白进行聚类 (默认30%序列一致性)
2. 构建共现矩阵: 行=RT, 列=蛋白簇, 值=0/1 (存在/不存在)
3. 计算每个蛋白簇的出现频率统计

与DefenseFinder方法的区别:
- DefenseFinder: 只检测已知的防御系统
- MMseqs2聚类: 无监督聚类，可以发现未知的关联蛋白

输入: Step6的all_neighbors.faa和neighbor_matrix.tsv
输出:
  - clusters_rep.faa: 聚类代表序列
  - cluster_members.tsv: 聚类成员关系
  - rt_cluster_matrix.tsv: RT-蛋白簇共现矩阵 (0/1)
  - cluster_stats.tsv: 蛋白簇统计信息
"""

import argparse
import pandas as pd
import subprocess
import sys
import logging
import shutil
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


def check_mmseqs():
    """检查MMseqs2是否可用"""
    return shutil.which('mmseqs') is not None


def run_mmseqs_cluster(input_fasta, output_prefix, tmp_dir, min_seq_id=0.3, coverage=0.8,
                       cov_mode=0, threads=8, logger=None):
    """
    运行MMseqs2聚类

    参数:
        min_seq_id: 最小序列一致性 (默认0.3 = 30%)
        coverage: 最小覆盖度 (默认0.8 = 80%)
        cov_mode: 覆盖度模式 (0=双向, 1=目标, 2=查询)
    """
    cmd = [
        "mmseqs", "easy-cluster",
        str(input_fasta),
        str(output_prefix),
        str(tmp_dir),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "--threads", str(threads),
        "-v", "2"  # 适中的详细程度
    ]

    logger.info(f"运行MMseqs2聚类...")
    logger.debug(f"命令: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200  # 2小时超时
        )

        if result.returncode != 0:
            logger.error(f"MMseqs2错误: {result.stderr}")
            return False

        logger.info("MMseqs2聚类完成")
        return True

    except subprocess.TimeoutExpired:
        logger.error("MMseqs2超时")
        return False
    except Exception as e:
        logger.error(f"MMseqs2运行失败: {e}")
        return False


def parse_mmseqs_clusters(cluster_tsv, logger):
    """
    解析MMseqs2聚类结果

    MMseqs2 easy-cluster 输出格式 (TSV):
    - 列1: 代表序列ID (cluster representative)
    - 列2: 成员序列ID (cluster member)

    返回:
        cluster_map: {member_id: representative_id}
        cluster_members: {representative_id: [member_ids]}
    """
    cluster_map = {}  # member -> representative
    cluster_members = defaultdict(list)  # representative -> [members]

    with open(cluster_tsv, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) >= 2:
                rep_id = parts[0]
                member_id = parts[1]

                cluster_map[member_id] = rep_id
                cluster_members[rep_id].append(member_id)

    logger.info(f"解析到 {len(cluster_members)} 个蛋白簇")
    logger.info(f"总成员数: {len(cluster_map)}")

    # 统计簇大小分布
    sizes = [len(members) for members in cluster_members.values()]
    logger.info(f"簇大小范围: {min(sizes)} - {max(sizes)}")
    logger.info(f"单成员簇: {sum(1 for s in sizes if s == 1)}")
    logger.info(f"多成员簇 (>1): {sum(1 for s in sizes if s > 1)}")
    logger.info(f"大簇 (>5): {sum(1 for s in sizes if s > 5)}")

    return cluster_map, dict(cluster_members)


def extract_rt_from_neighbor_id(neighbor_id):
    """
    从邻近基因ID中提取RT ID

    邻近基因ID格式: RT_ID|locus_tag|position
    例如: genome_contig_123-456|ctg1_001|upstream_5000bp
    """
    parts = neighbor_id.split('|')
    if len(parts) >= 1:
        return parts[0]
    return neighbor_id


def build_cooccurrence_matrix(neighbor_matrix_file, cluster_map, cluster_members, logger):
    """
    构建RT-蛋白簇共现矩阵

    矩阵结构:
    - 行: 每个RT (rt_id)
    - 列: 每个蛋白簇 (representative_id)
    - 值: 1 如果该RT的邻近区域包含该簇的成员, 否则为0

    返回:
        matrix_df: 共现矩阵 (DataFrame)
        rt_cluster_counts: 每个RT-簇组合的成员计数
    """
    # 读取邻近基因信息
    neighbor_df = pd.read_csv(neighbor_matrix_file, sep='\t')
    logger.info(f"读取 {len(neighbor_df)} 条邻近基因记录")

    # 提取所有RT ID
    rt_ids = neighbor_df['rt_id'].unique().tolist()
    logger.info(f"RT数量: {len(rt_ids)}")

    # 提取所有簇代表ID
    cluster_reps = list(cluster_members.keys())
    logger.info(f"蛋白簇数量: {len(cluster_reps)}")

    # 构建RT -> 邻近基因ID的映射
    rt_to_neighbors = defaultdict(set)
    for _, row in neighbor_df.iterrows():
        rt_id = row['rt_id']
        neighbor_id = row['neighbor_id']
        rt_to_neighbors[rt_id].add(neighbor_id)

    # 构建共现矩阵
    logger.info("构建共现矩阵...")

    # 存储矩阵数据
    matrix_data = {}
    rt_cluster_counts = []  # 详细计数信息

    for rt_id in rt_ids:
        neighbor_ids = rt_to_neighbors[rt_id]
        matrix_data[rt_id] = {}

        for cluster_rep in cluster_reps:
            # 检查该RT的邻近基因是否属于这个簇
            cluster_member_set = set(cluster_members[cluster_rep])
            intersection = neighbor_ids & cluster_member_set

            if intersection:
                matrix_data[rt_id][cluster_rep] = 1
                # 记录详细信息
                rt_cluster_counts.append({
                    'rt_id': rt_id,
                    'cluster_rep': cluster_rep,
                    'member_count': len(intersection),
                    'members': ';'.join(list(intersection)[:5])  # 只记录前5个
                })
            else:
                matrix_data[rt_id][cluster_rep] = 0

    # 转换为DataFrame
    matrix_df = pd.DataFrame.from_dict(matrix_data, orient='index')
    matrix_df.index.name = 'rt_id'

    # 计算统计信息
    cluster_occurrence = matrix_df.sum(axis=0)  # 每个簇在多少个RT中出现
    rt_diversity = matrix_df.sum(axis=1)  # 每个RT有多少个不同的簇

    logger.info(f"矩阵维度: {matrix_df.shape[0]} RTs × {matrix_df.shape[1]} 蛋白簇")
    logger.info(f"非零条目: {(matrix_df > 0).sum().sum()}")
    logger.info(f"矩阵稀疏度: {1 - (matrix_df > 0).sum().sum() / (matrix_df.shape[0] * matrix_df.shape[1]):.2%}")

    return matrix_df, pd.DataFrame(rt_cluster_counts)


def calculate_cluster_stats(matrix_df, cluster_members, min_occurrence=1, logger=None):
    """
    计算蛋白簇统计信息

    统计项:
    - occurrence_count: 在多少个RT中出现
    - occurrence_rate: 出现率 (占总RT数的比例)
    - cluster_size: 簇成员数量
    """
    total_rts = len(matrix_df)
    cluster_stats = []

    for cluster_rep in matrix_df.columns:
        occurrence_count = matrix_df[cluster_rep].sum()
        occurrence_rate = occurrence_count / total_rts

        cluster_stats.append({
            'cluster_rep': cluster_rep,
            'cluster_size': len(cluster_members.get(cluster_rep, [])),
            'occurrence_count': int(occurrence_count),
            'occurrence_rate': round(occurrence_rate, 4),
            'in_rt_count': int(occurrence_count),  # 同义字段，便于理解
        })

    stats_df = pd.DataFrame(cluster_stats)
    stats_df = stats_df.sort_values('occurrence_count', ascending=False)

    # 筛选出现次数大于阈值的簇
    significant_clusters = stats_df[stats_df['occurrence_count'] >= min_occurrence]

    if logger:
        logger.info(f"\n蛋白簇统计:")
        logger.info(f"  总蛋白簇数: {len(stats_df)}")
        logger.info(f"  出现次数 >= {min_occurrence} 的簇: {len(significant_clusters)}")
        if len(significant_clusters) > 0:
            logger.info(f"  最高出现次数: {significant_clusters['occurrence_count'].max()}")
            logger.info(f"  最高出现率: {significant_clusters['occurrence_rate'].max():.1%}")

    return stats_df


def main():
    parser = argparse.ArgumentParser(
        description="Step 7: MMseqs2聚类与共现矩阵构建 (Mestre方法)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法
  python step07_mmseqs_cluster.py -i 06_neighbors -o 07_clustering

  # 自定义聚类参数
  python step07_mmseqs_cluster.py -i 06_neighbors -o 07_clustering \\
    --min-seq-id 0.3 --coverage 0.8

  # 使用已有的聚类结果
  python step07_mmseqs_cluster.py -i 06_neighbors -o 07_clustering \\
    --existing-clusters existing_cluster.tsv

方法说明:
  本步骤采用Mestre et al., 2020的方法:
  1. 使用MMseqs2对所有邻近蛋白进行无监督聚类
  2. 构建RT-蛋白簇的共现矩阵 (0/1二值矩阵)
  3. 统计每个蛋白簇的出现频率

  这为后续的Phyvalue统计关联分析奠定基础。
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Step6输出目录 (包含all_neighbors.faa和neighbor_matrix.tsv)')
    parser.add_argument('-o', '--output', default='07_clustering',
                        help='输出目录 (默认: 07_clustering)')

    # MMseqs2参数
    parser.add_argument('--min-seq-id', type=float, default=0.3,
                        help='最小序列一致性 (默认: 0.3 = 30%%)')
    parser.add_argument('--coverage', type=float, default=0.8,
                        help='最小覆盖度 (默认: 0.8 = 80%%)')
    parser.add_argument('--cov-mode', type=int, default=0, choices=[0, 1, 2],
                        help='覆盖度模式: 0=双向, 1=目标, 2=查询 (默认: 0)')
    parser.add_argument('--threads', type=int, default=8,
                        help='线程数 (默认: 8)')

    # 可选：使用已有的聚类结果
    parser.add_argument('--existing-clusters',
                        help='使用已有的MMseqs2聚类结果TSV文件')

    parser.add_argument('--min-occurrence', type=int, default=1,
                        help='统计时的最小出现次数阈值 (默认: 1)')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    # 检查MMseqs2
    if not args.existing_clusters and not check_mmseqs():
        logger.error("MMseqs2未安装或不在PATH中")
        logger.error("安装方法: conda install -c bioconda mmseqs2")
        sys.exit(1)

    # 验证输入
    input_dir = Path(args.input)
    neighbors_fasta = input_dir / "all_neighbors.faa"
    neighbor_matrix = input_dir / "neighbor_matrix.tsv"

    if not neighbors_fasta.exists():
        logger.error(f"邻近蛋白文件不存在: {neighbors_fasta}")
        sys.exit(1)
    if not neighbor_matrix.exists():
        logger.error(f"邻近基因矩阵不存在: {neighbor_matrix}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = output_dir / "tmp_mmseqs"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 7: MMseqs2聚类与共现矩阵构建 (Mestre方法)")
    logger.info("=" * 60)
    logger.info(f"输入目录: {input_dir}")
    logger.info(f"输出目录: {output_dir}")

    # 运行MMseqs2聚类或使用已有结果
    cluster_prefix = output_dir / "clusters"

    if args.existing_clusters:
        cluster_tsv = Path(args.existing_clusters)
        if not cluster_tsv.exists():
            logger.error(f"聚类结果文件不存在: {cluster_tsv}")
            sys.exit(1)
        logger.info(f"使用已有聚类结果: {cluster_tsv}")
    else:
        logger.info(f"\n聚类参数:")
        logger.info(f"  序列一致性: {args.min_seq_id:.0%}")
        logger.info(f"  覆盖度: {args.coverage:.0%}")
        logger.info(f"  覆盖度模式: {args.cov_mode}")
        logger.info(f"  线程数: {args.threads}")

        success = run_mmseqs_cluster(
            neighbors_fasta,
            cluster_prefix,
            tmp_dir,
            min_seq_id=args.min_seq_id,
            coverage=args.coverage,
            cov_mode=args.cov_mode,
            threads=args.threads,
            logger=logger
        )

        if not success:
            logger.error("MMseqs2聚类失败")
            sys.exit(1)

        cluster_tsv = Path(f"{cluster_prefix}_cluster.tsv")

    # 解析聚类结果
    logger.info("\n解析聚类结果...")
    cluster_map, cluster_members = parse_mmseqs_clusters(cluster_tsv, logger)

    # 保存聚类成员关系
    cluster_member_list = []
    for rep, members in cluster_members.items():
        for member in members:
            cluster_member_list.append({
                'cluster_rep': rep,
                'member_id': member,
                'is_representative': member == rep
            })

    member_df = pd.DataFrame(cluster_member_list)
    member_file = output_dir / "cluster_members.tsv"
    member_df.to_csv(member_file, sep='\t', index=False)
    logger.info(f"✓ 聚类成员关系: {member_file}")

    # 构建共现矩阵
    logger.info("\n构建RT-蛋白簇共现矩阵...")
    matrix_df, rt_cluster_counts_df = build_cooccurrence_matrix(
        neighbor_matrix, cluster_map, cluster_members, logger
    )

    # 保存共现矩阵
    matrix_file = output_dir / "rt_cluster_matrix.tsv"
    matrix_df.to_csv(matrix_file, sep='\t')
    logger.info(f"✓ 共现矩阵: {matrix_file}")

    # 保存详细计数
    if not rt_cluster_counts_df.empty:
        counts_file = output_dir / "rt_cluster_counts.tsv"
        rt_cluster_counts_df.to_csv(counts_file, sep='\t', index=False)
        logger.info(f"✓ RT-簇详细计数: {counts_file}")

    # 计算蛋白簇统计
    logger.info("\n计算蛋白簇统计...")
    stats_df = calculate_cluster_stats(
        matrix_df, cluster_members,
        min_occurrence=args.min_occurrence,
        logger=logger
    )

    stats_file = output_dir / "cluster_stats.tsv"
    stats_df.to_csv(stats_file, sep='\t', index=False)
    logger.info(f"✓ 蛋白簇统计: {stats_file}")

    # 输出高频簇信息
    high_freq_clusters = stats_df[stats_df['occurrence_count'] >= 5]
    if not high_freq_clusters.empty:
        high_freq_file = output_dir / "high_frequency_clusters.tsv"
        high_freq_clusters.to_csv(high_freq_file, sep='\t', index=False)
        logger.info(f"✓ 高频簇 (>=5 RTs): {high_freq_file} ({len(high_freq_clusters)} 个)")

    # 清理临时目录
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir, ignore_errors=True)

    # 总结
    logger.info("\n" + "=" * 60)
    logger.info("聚类与矩阵构建完成!")
    logger.info(f"  RT数量: {matrix_df.shape[0]}")
    logger.info(f"  蛋白簇数量: {matrix_df.shape[1]}")
    logger.info(f"  高频簇 (>=5 RTs): {len(high_freq_clusters)}")
    logger.info(f"  输出目录: {output_dir}")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
