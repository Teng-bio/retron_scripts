#!/usr/bin/env python3
"""
Step 3: 排除非Retron RT

基于多种特征排除其他类型的RT：
1. Group II Intron RT - 含有内切酶结构域，通常位于基因内部
2. DGR RT - 与Variable Repeat (VR)和Template Repeat (TR)相关
3. CRISPR-associated RT - 与Cas基因连锁

检测方法：
- 邻近基因产物关键词匹配
- 结构域特征模式匹配
- 基因组位置特征

参考: Zimmerly & Wu, 2015; Millman et al., 2020
"""

import argparse
import pandas as pd
import re
import sys
import logging
from pathlib import Path


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


# 排除关键词列表
GROUP_II_KEYWORDS = [
    'group II intron', 'intron maturase', 'maturase',
    'ltrA', 'splicing', 'self-splicing',
    'intron-encoded', 'retroelement'
]

DGR_KEYWORDS = [
    'diversity-generating', 'DGR', 'variable repeat',
    'template repeat', 'avd', 'accessory variability determinant',
    'target gene', 'MTase'  # DGR常伴随的甲基转移酶
]

CRISPR_KEYWORDS = [
    'cas1', 'cas2', 'cas3', 'cas9', 'cas12', 'cas13',
    'crispr', 'spacer', 'repeat', 'adaptation'
]

# 可能是Retron效应蛋白的关键词（正面证据）
EFFECTOR_KEYWORDS = [
    'nuclease', 'endonuclease', 'HNH', 'OLD',
    'ATPase', 'transmembrane', 'toxin', 'antitoxin',
    'ribosyltransferase', 'cold shock'
]


def check_group_ii_intron(product_description, sequence=None):
    """
    检查是否为Group II Intron RT

    特征：
    - 产物描述含有intron/maturase关键词
    - RT通常较长(>500aa)且含有内切酶结构域
    """
    # 处理空值和NaN
    if not product_description or (isinstance(product_description, float) and pd.isna(product_description)):
        return False, []

    product_lower = str(product_description).lower()
    matches = []

    for kw in GROUP_II_KEYWORDS:
        if kw.lower() in product_lower:
            matches.append(kw)

    return len(matches) > 0, matches


def check_dgr(product_description, neighbor_products=None):
    """
    检查是否为DGR RT

    特征：
    - 产物描述含有DGR相关关键词
    - 邻近基因含有variable repeat或template repeat
    """
    matches = []

    # 处理空值和NaN
    if product_description and not (isinstance(product_description, float) and pd.isna(product_description)):
        product_lower = str(product_description).lower()
        for kw in DGR_KEYWORDS:
            if kw.lower() in product_lower:
                matches.append(f"self:{kw}")

    if neighbor_products:
        for np in neighbor_products:
            # 处理空值和NaN
            if not np or (isinstance(np, float) and pd.isna(np)):
                continue
            np_lower = str(np).lower()
            for kw in DGR_KEYWORDS:
                if kw.lower() in np_lower:
                    matches.append(f"neighbor:{kw}")

    return len(matches) > 0, matches


def check_crispr_associated(product_description, neighbor_products=None):
    """
    检查是否为CRISPR-associated RT

    特征：
    - 邻近基因含有Cas蛋白
    - 产物描述含有CRISPR相关关键词
    """
    matches = []

    # 处理空值和NaN
    if product_description and not (isinstance(product_description, float) and pd.isna(product_description)):
        product_lower = str(product_description).lower()
        for kw in CRISPR_KEYWORDS:
            if kw.lower() in product_lower:
                matches.append(f"self:{kw}")

    if neighbor_products:
        for np in neighbor_products:
            # 处理空值和NaN
            if not np or (isinstance(np, float) and pd.isna(np)):
                continue
            np_lower = str(np).lower()
            for kw in CRISPR_KEYWORDS:
                if kw.lower() in np_lower:
                    matches.append(f"neighbor:{kw}")

    return len(matches) > 0, matches


def check_effector_nearby(neighbor_products):
    """
    检查邻近是否有效应蛋白（正面证据）
    """
    if not neighbor_products:
        return False, []

    matches = []
    for np in neighbor_products:
        # 处理空值和NaN
        if not np or (isinstance(np, float) and pd.isna(np)):
            continue
        np_lower = str(np).lower()
        for kw in EFFECTOR_KEYWORDS:
            if kw.lower() in np_lower:
                matches.append(f"{kw}")

    return len(matches) > 0, list(set(matches))


def classify_rt(row, neighbor_info=None):
    """
    对单个RT进行分类

    返回: (classification, exclusion_reason, evidence)
    """
    product = row.get('antismash_product', '') or row.get('product', '') or ''
    neighbor_products = []

    if neighbor_info:
        neighbor_products = neighbor_info

    # 检查各种非Retron类型
    is_group_ii, group_ii_evidence = check_group_ii_intron(product)
    is_dgr, dgr_evidence = check_dgr(product, neighbor_products)
    is_crispr, crispr_evidence = check_crispr_associated(product, neighbor_products)
    has_effector, effector_evidence = check_effector_nearby(neighbor_products)

    # 判定逻辑
    if is_group_ii:
        return 'group_ii_intron', 'Group II intron特征', group_ii_evidence
    elif is_crispr:
        return 'crispr_associated', 'CRISPR-associated RT', crispr_evidence
    elif is_dgr:
        return 'dgr', 'DGR RT', dgr_evidence
    elif has_effector:
        return 'retron_candidate', 'Retron候选(有效应蛋白证据)', effector_evidence
    else:
        return 'unclassified', '无明确排除证据', []


def main():
    parser = argparse.ArgumentParser(
        description="Step 3: 排除非Retron RT (Group II/DGR/CRISPR)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Step2的motif过滤结果')
    parser.add_argument('-n', '--neighbors',
                        help='邻近基因信息文件 (可选，用于更准确分类)')
    parser.add_argument('-o', '--output', default='03_classified',
                        help='输出目录')
    parser.add_argument('--keep-excluded', action='store_true',
                        help='保留被排除的RT（用于审查）')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"输入文件不存在: {input_file}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 3: RT分类与非Retron排除")
    logger.info("=" * 60)

    # 读取数据
    df = pd.read_csv(input_file, sep='\t')
    logger.info(f"读取到 {len(df)} 条RT候选")

    # 加载邻近基因信息（如果有）
    neighbor_info = {}
    if args.neighbors and Path(args.neighbors).exists():
        neighbors_df = pd.read_csv(args.neighbors, sep='\t')
        for _, row in neighbors_df.iterrows():
            key = f"{row['genome']}|{row.get('query', '')}"
            if key not in neighbor_info:
                neighbor_info[key] = []
            neighbor_info[key].append(row.get('product', ''))
        logger.info(f"加载了 {len(neighbor_info)} 个候选的邻近基因信息")

    # 分类每个RT
    logger.info("\n分类RT...")
    results = []

    for _, row in df.iterrows():
        # 获取邻近基因信息
        key = f"{row['genome']}|{row.get('query', '')}"
        neighbors = neighbor_info.get(key, [])

        classification, reason, evidence = classify_rt(row, neighbors)

        result = row.to_dict()
        result.update({
            'rt_classification': classification,
            'classification_reason': reason,
            'classification_evidence': ';'.join(evidence) if evidence else ''
        })
        results.append(result)

    result_df = pd.DataFrame(results)

    # 统计
    stats = result_df['rt_classification'].value_counts()
    logger.info("\nRT分类统计:")
    for cls, count in stats.items():
        logger.info(f"  {cls}: {count}")

    # 分离Retron候选和排除的
    excluded_types = ['group_ii_intron', 'crispr_associated', 'dgr']
    retron_df = result_df[~result_df['rt_classification'].isin(excluded_types)]
    excluded_df = result_df[result_df['rt_classification'].isin(excluded_types)]

    logger.info(f"\nRetron候选: {len(retron_df)} 条")
    logger.info(f"排除: {len(excluded_df)} 条")

    # 保存结果
    # Retron候选
    retron_output = output_dir / "retron_candidates.tsv"
    retron_df.to_csv(retron_output, sep='\t', index=False)
    logger.info(f"\n✓ Retron候选: {retron_output}")

    # 完整分类结果
    full_output = output_dir / "classification_full.tsv"
    result_df.to_csv(full_output, sep='\t', index=False)
    logger.info(f"✓ 完整分类: {full_output}")

    # 排除的（如果需要）
    if args.keep_excluded and not excluded_df.empty:
        excluded_output = output_dir / "excluded_rt.tsv"
        excluded_df.to_csv(excluded_output, sep='\t', index=False)
        logger.info(f"✓ 排除的RT: {excluded_output}")

    # 统计报告
    stats_file = output_dir / "classification_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("RT分类统计\n")
        f.write("=" * 40 + "\n")
        f.write(f"输入: {len(df)} 条\n\n")
        f.write("分类结果:\n")
        for cls, count in stats.items():
            pct = count / len(df) * 100
            f.write(f"  {cls}: {count} ({pct:.1f}%)\n")
        f.write(f"\n保留为Retron候选: {len(retron_df)} 条\n")
        f.write(f"排除: {len(excluded_df)} 条\n")

    logger.info("\n" + "=" * 60)
    logger.info(f"分类完成! 保留 {len(retron_df)}/{len(df)} 条")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
