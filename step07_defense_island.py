#!/usr/bin/env python3
"""
Step 7: 防御岛分析

使用DefenseFinder检测基因簇邻近基因中的防御系统。
分析Retron候选是否位于防御岛内（即邻近基因中是否存在防御系统）。

约60%的Retron位于已知防御系统附近（Millman et al., 2020）

支持两种模式:
1. 直接运行DefenseFinder（需要在DefenseFinder环境中）
2. 使用已有的DefenseFinder结果（通过 --df-results 参数）

如果DefenseFinder需要在单独环境中运行，请先使用 run_defensefinder_batch.py
"""

import argparse
import pandas as pd
import subprocess
import sys
import logging
import shutil
import tempfile
from pathlib import Path


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


def check_defensefinder():
    """检查DefenseFinder是否可用"""
    return shutil.which('defense-finder') is not None


def run_defensefinder_on_cluster(fasta_file, output_dir, logger):
    """
    对单个基因簇的邻近蛋白运行DefenseFinder

    返回: (systems_df, genes_df) 或 (None, None)
    """
    try:
        cmd = [
            'defense-finder', 'run',
            str(fasta_file),
            '-o', str(output_dir),
            '--db-type', 'unordered'
        ]

        logger.debug(f"运行命令: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )

        if result.returncode != 0 and 'error' in result.stderr.lower():
            logger.debug(f"DefenseFinder警告: {result.stderr[:200]}")

        # 查找结果文件
        systems_file = output_dir / "defense_finder_systems.tsv"
        genes_file = output_dir / "defense_finder_genes.tsv"

        systems_df = None
        genes_df = None

        if systems_file.exists():
            try:
                systems_df = pd.read_csv(systems_file, sep='\t')
            except:
                pass
        if genes_file.exists():
            try:
                genes_df = pd.read_csv(genes_file, sep='\t')
            except:
                pass

        return systems_df, genes_df

    except subprocess.TimeoutExpired:
        logger.debug("DefenseFinder超时")
        return None, None
    except Exception as e:
        logger.debug(f"DefenseFinder运行失败: {e}")
        return None, None


def parse_existing_results(results_dir, cluster_id, logger):
    """解析已有的DefenseFinder结果"""
    results_dir = Path(results_dir)
    cluster_results = results_dir / cluster_id

    genes_file = cluster_results / "defense_finder_genes.tsv"

    if not genes_file.exists():
        return None, None

    try:
        genes_df = pd.read_csv(genes_file, sep='\t')
        systems_file = cluster_results / "defense_finder_systems.tsv"
        systems_df = None
        if systems_file.exists():
            systems_df = pd.read_csv(systems_file, sep='\t')
        return systems_df, genes_df
    except Exception as e:
        logger.debug(f"解析结果失败 {cluster_id}: {e}")
        return None, None


def analyze_clusters(clusters_dir, cluster_summary_file, output_dir, logger, df_results_dir=None):
    """
    分析所有基因簇的防御系统

    参数:
        df_results_dir: 已有的DefenseFinder结果目录（可选）

    返回: 带有防御系统信息的DataFrame
    """
    clusters_dir = Path(clusters_dir)
    output_dir = Path(output_dir)

    # 读取基因簇汇总
    summary_df = pd.read_csv(cluster_summary_file, sep='\t')
    logger.info(f"读取到 {len(summary_df)} 个基因簇")

    use_existing = df_results_dir is not None
    if use_existing:
        logger.info(f"使用已有DefenseFinder结果: {df_results_dir}")

    results = []
    defense_details = []

    for idx, row in summary_df.iterrows():
        cluster_id = row['cluster_id']
        cluster_path = clusters_dir / cluster_id

        # 进度报告
        if (idx + 1) % 50 == 0 or idx == 0:
            logger.info(f"处理进度: {idx + 1}/{len(summary_df)}")

        # 获取DefenseFinder结果
        systems_df, genes_df = None, None

        if use_existing:
            # 使用已有结果
            systems_df, genes_df = parse_existing_results(df_results_dir, cluster_id, logger)
        else:
            # 运行DefenseFinder
            fasta_file = cluster_path / "neighbor_proteins.fasta"
            if fasta_file.exists():
                with tempfile.TemporaryDirectory() as tmp_dir:
                    tmp_path = Path(tmp_dir)
                    systems_df, genes_df = run_defensefinder_on_cluster(fasta_file, tmp_path, logger)

        # 解析结果
        has_defense = False
        defense_systems = []
        defense_genes = []

        if genes_df is not None and not genes_df.empty:
            has_defense = True
            # 提取防御系统信息
            if 'type' in genes_df.columns:
                defense_systems = genes_df['type'].unique().tolist()
            elif 'sys_id' in genes_df.columns:
                defense_systems = genes_df['sys_id'].unique().tolist()

            if 'hit_id' in genes_df.columns:
                defense_genes = genes_df['hit_id'].tolist()
            elif 'gene_name' in genes_df.columns:
                defense_genes = genes_df['gene_name'].tolist()

            # 记录详细信息
            for _, gene_row in genes_df.iterrows():
                defense_details.append({
                    'cluster_id': cluster_id,
                    'genome': row['genome'],
                    'contig': row['contig'],
                    **gene_row.to_dict()
                })

            logger.debug(f"  {cluster_id}: 发现 {len(defense_genes)} 个防御基因")

        results.append({
            **row.to_dict(),
            'has_defense_system': has_defense,
            'defense_systems': ';'.join(str(s) for s in defense_systems),
            'defense_genes': ';'.join(str(g) for g in defense_genes),
            'defense_gene_count': len(defense_genes)
        })

    return pd.DataFrame(results), pd.DataFrame(defense_details) if defense_details else None


def main():
    parser = argparse.ArgumentParser(
        description="Step 7: 防御岛分析 (DefenseFinder)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 模式1: 直接运行DefenseFinder（需要在DefenseFinder环境中）
  python step07_defense_island.py -i 06_neighbors -o 07_defense_island

  # 模式2: 使用已有的DefenseFinder结果
  # 先在DefenseFinder环境中运行:
  #   python run_defensefinder_batch.py -i 06_neighbors/clusters -o df_results
  # 然后:
  python step07_defense_island.py -i 06_neighbors --df-results df_results -o 07_defense_island

注意:
  - 运行前请先执行: defense-finder update
  - 输入应为step06生成的输出目录（包含clusters子目录）
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Step06输出目录 (包含cluster_summary.tsv和clusters/)')
    parser.add_argument('-o', '--output', default='07_defense_island',
                        help='输出目录')
    parser.add_argument('--df-results',
                        help='已有的DefenseFinder结果目录 (由run_defensefinder_batch.py生成)')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    # 检查DefenseFinder（仅在不使用已有结果时）
    if not args.df_results and not check_defensefinder():
        logger.error("DefenseFinder未安装或不在PATH中")
        logger.error("选项1: 安装DefenseFinder: conda install -c bioconda defense-finder")
        logger.error("选项2: 在DefenseFinder环境中运行 run_defensefinder_batch.py，然后使用 --df-results")
        sys.exit(1)

    input_dir = Path(args.input)

    # 支持两种输入格式
    if input_dir.is_file():
        parent_dir = input_dir.parent.parent / "06_neighbors"
        if parent_dir.exists():
            input_dir = parent_dir
        else:
            logger.error("输入应为step06的输出目录，包含cluster_summary.tsv和clusters/")
            sys.exit(1)

    cluster_summary_file = input_dir / "cluster_summary.tsv"
    clusters_dir = input_dir / "clusters"

    if not cluster_summary_file.exists():
        logger.error(f"未找到基因簇汇总文件: {cluster_summary_file}")
        logger.error("请先运行step06生成基因簇数据")
        sys.exit(1)

    if not clusters_dir.exists():
        logger.error(f"未找到基因簇目录: {clusters_dir}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 7: 防御岛分析 (基于邻近基因)")
    logger.info("=" * 60)
    logger.info(f"基因簇目录: {clusters_dir}")
    if args.df_results:
        logger.info(f"使用已有结果: {args.df_results}")

    # 分析防御系统
    df_results_dir = Path(args.df_results) if args.df_results else None
    result_df, details_df = analyze_clusters(
        clusters_dir, cluster_summary_file, output_dir, logger, df_results_dir
    )

    # 统计
    has_defense_count = result_df['has_defense_system'].sum()
    total_count = len(result_df)
    logger.info(f"\n含防御系统的基因簇: {has_defense_count}/{total_count} "
                f"({has_defense_count/total_count*100:.1f}%)")

    # 保存结果
    output_file = output_dir / "defense_analysis.tsv"
    result_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"\n✓ 结果: {output_file}")

    # 保存防御基因详细信息
    if details_df is not None and not details_df.empty:
        details_file = output_dir / "defense_genes_detail.tsv"
        details_df.to_csv(details_file, sep='\t', index=False)
        logger.info(f"✓ 防御基因详情: {details_file}")

    # 含防御系统的基因簇
    defense_clusters = result_df[result_df['has_defense_system'] == True]
    if not defense_clusters.empty:
        defense_file = output_dir / "clusters_with_defense.tsv"
        defense_clusters.to_csv(defense_file, sep='\t', index=False)
        logger.info(f"✓ 含防御系统基因簇: {defense_file} ({len(defense_clusters)} 个)")

    # 统计报告
    stats_file = output_dir / "defense_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("防御岛分析统计\n")
        f.write("=" * 40 + "\n")
        f.write(f"分析基因簇数: {total_count}\n")
        f.write(f"含防御系统: {has_defense_count} ({has_defense_count/total_count*100:.1f}%)\n")

        if 'defense_systems' in result_df.columns:
            all_systems = []
            for s in result_df['defense_systems'].dropna():
                if s:
                    all_systems.extend(s.split(';'))
            if all_systems:
                system_counts = pd.Series(all_systems).value_counts()
                f.write(f"\n发现的防御系统类型:\n")
                for sys_name, count in system_counts.items():
                    f.write(f"  {sys_name}: {count}\n")

    logger.info("\n" + "=" * 60)
    logger.info("分析完成!")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
