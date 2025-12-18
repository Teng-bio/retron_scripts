#!/usr/bin/env python3
"""
Step 9: 综合报告生成

汇总所有分析结果，生成最终的Retron候选报告：
- 整合各步骤结果
- 计算综合置信度评分
- 生成分级候选列表
"""

import argparse
import pandas as pd
import sys
import logging
from pathlib import Path
from datetime import datetime


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


def calculate_confidence_score(row):
    """
    计算综合置信度评分 (0-100)

    评分维度：
    - Motif (30分): NAXXH + VTG
    - RNA结构 (25分): 茎长度、分支G、MFE
    - 效应蛋白 (25分): 邻近效应蛋白
    - 防御岛 (20分): 位于防御系统附近
    """
    score = 0
    evidence = []

    # 1. Motif评分 (30分)
    motif_class = row.get('motif_classification', '')
    if motif_class == 'high_confidence':
        score += 30
        evidence.append("Motif:NAXXH+VTG")
    elif motif_class == 'medium_confidence':
        score += 20
        evidence.append("Motif:NAXXH")
    elif motif_class == 'low_confidence':
        score += 10
        evidence.append("Motif:VTG_only")

    # 2. RNA结构评分 (25分)
    # 茎长度
    stem = row.get('longest_stem', 0)
    if pd.notna(stem):
        stem = int(stem)
        if stem >= 20:
            score += 15
            evidence.append(f"Stem:{stem}bp")
        elif stem >= 15:
            score += 10
            evidence.append(f"Stem:{stem}bp")
        elif stem >= 10:
            score += 5

    # 分支G
    branch_g = row.get('branch_g_count', 0)
    if pd.notna(branch_g) and int(branch_g) > 0:
        score += 5
        evidence.append(f"BranchG:{int(branch_g)}")

    # MFE
    mfe = row.get('mfe', 0)
    if pd.notna(mfe):
        mfe = float(mfe)
        if mfe < -25:
            score += 5
        elif mfe < -20:
            score += 3

    # 3. 效应蛋白评分 (25分)
    has_effector = row.get('has_effector', False)
    if has_effector:
        score += 15
        eff_types = row.get('effector_types', '')
        if eff_types:
            evidence.append(f"Effector:{eff_types.split(';')[0]}")

        # 额外加分：特定类型效应蛋白
        if any(t in str(eff_types).lower() for t in ['hnh', 'nuclease', 'atpase', '2tm']):
            score += 10

    # 4. 防御岛评分 (20分)
    in_defense = row.get('in_defense_island', False)
    if in_defense:
        score += 20
        evidence.append("DefenseIsland")

    return min(score, 100), ';'.join(evidence)


def classify_confidence(score):
    """根据分数分类置信度"""
    if score >= 70:
        return 'HIGH'
    elif score >= 45:
        return 'MEDIUM'
    elif score >= 25:
        return 'LOW'
    else:
        return 'VERY_LOW'


def merge_results(files, logger):
    """
    合并各步骤的结果

    核心逻辑：
    1. 以classified结果为RT蛋白主表（每个RT蛋白一条记录）
    2. 从RNA结果中为每个RT选择最佳窗口（MFE最低）
    3. 合并motif、defense、effector信息

    合并键：genome + protein_id（RT蛋白唯一标识）
    """
    dfs = {}

    for name, path in files.items():
        if path and Path(path).exists():
            try:
                dfs[name] = pd.read_csv(path, sep='\t')
                logger.info(f"  加载 {name}: {len(dfs[name])} 条")
            except Exception as e:
                logger.warning(f"  加载 {name} 失败: {e}")

    if not dfs:
        logger.error("没有可用的结果文件")
        return None

    # 以classified为RT蛋白主表
    if 'classified' not in dfs:
        logger.error("缺少classified结果，无法确定RT蛋白列表")
        return None

    # classified中的contig列实际是protein_id
    rt_df = dfs['classified'].copy()
    rt_df = rt_df.rename(columns={'contig': 'protein_id'})
    logger.info(f"\nRT蛋白主表: {len(rt_df)} 个RT蛋白")

    # 处理RNA结果：为每个RT选择最佳窗口
    if 'rna' in dfs:
        rna_df = dfs['rna'].copy()

        # 检查是否有protein_id列
        if 'protein_id' in rna_df.columns:
            # 按genome+protein_id分组，选择MFE最低的窗口
            rna_df['mfe'] = pd.to_numeric(rna_df['mfe'], errors='coerce')
            best_rna = rna_df.loc[rna_df.groupby(['genome', 'protein_id'])['mfe'].idxmin()]
            logger.info(f"  RNA窗口去重: {len(rna_df)} → {len(best_rna)} (每个RT选最佳)")
        else:
            # 旧格式：尝试从seq_id解析protein_id
            def extract_protein_id(seq_id):
                parts = str(seq_id).split('|')
                if len(parts) >= 6:  # 新格式
                    return parts[2]
                return ''

            rna_df['protein_id'] = rna_df['seq_id'].apply(extract_protein_id)
            if rna_df['protein_id'].str.len().sum() > 0:
                rna_df['mfe'] = pd.to_numeric(rna_df['mfe'], errors='coerce')
                best_rna = rna_df.loc[rna_df.groupby(['genome', 'protein_id'])['mfe'].idxmin()]
                logger.info(f"  RNA窗口去重: {len(rna_df)} → {len(best_rna)} (每个RT选最佳)")
            else:
                # 无法解析protein_id，按genome去重
                logger.warning("  RNA结果无protein_id，按genome去重（可能不准确）")
                best_rna = rna_df.loc[rna_df.groupby('genome')['mfe'].idxmin()]

        # 合并RNA信息到RT主表
        rna_cols = ['genome', 'protein_id', 'seq_id', 'sequence', 'structure', 'mfe',
                   'longest_stem', 'branch_g_count', 'branch_g_positions',
                   'global_start', 'global_end', 'strand']
        rna_cols = [c for c in rna_cols if c in best_rna.columns]

        result = rt_df.merge(
            best_rna[rna_cols],
            on=['genome', 'protein_id'],
            how='left',
            suffixes=('', '_rna')
        )
    else:
        result = rt_df

    # 合并motif信息（使用genome+protein_id）
    if 'motif' in dfs and 'motif_classification' not in result.columns:
        motif_df = dfs['motif'].copy()
        # motif中的contig也是protein_id
        if 'contig' in motif_df.columns:
            motif_df = motif_df.rename(columns={'contig': 'protein_id'})

        motif_cols = ['genome', 'protein_id', 'has_naxxh', 'has_vtg',
                      'motif_classification', 'naxxh_motifs', 'vtg_motifs']
        motif_cols = [c for c in motif_cols if c in motif_df.columns]

        if 'protein_id' in motif_cols:
            result = result.merge(
                motif_df[motif_cols].drop_duplicates(),
                on=['genome', 'protein_id'],
                how='left',
                suffixes=('', '_motif')
            )
        else:
            # 回退到只用genome
            motif_cols = [c for c in motif_cols if c != 'protein_id']
            result = result.merge(
                motif_df[motif_cols].drop_duplicates(),
                on=['genome'],
                how='left',
                suffixes=('', '_motif')
            )

    # 合并防御岛信息
    if 'defense' in dfs:
        defense_df = dfs['defense']
        if 'has_defense_system' in defense_df.columns:
            # 按genome合并（defense是基因簇级别的）
            defense_cols = ['genome', 'has_defense_system', 'defense_systems', 'defense_gene_count']
            defense_cols = [c for c in defense_cols if c in defense_df.columns]
            result = result.merge(
                defense_df[defense_cols].drop_duplicates(subset=['genome']),
                on=['genome'],
                how='left',
                suffixes=('', '_defense')
            )
            if 'has_defense_system' in result.columns and 'in_defense_island' not in result.columns:
                result['in_defense_island'] = result['has_defense_system']

    # 合并效应蛋白信息
    if 'effector' in dfs and 'has_effector' not in result.columns:
        effector_df = dfs['effector'].copy()

        # 检查是否有protein_id
        if 'protein_id' in effector_df.columns:
            effector_cols = ['genome', 'protein_id', 'has_effector', 'effector_count',
                           'effector_types', 'effector_details']
            effector_cols = [c for c in effector_cols if c in effector_df.columns]
            result = result.merge(
                effector_df[effector_cols].drop_duplicates(),
                on=['genome', 'protein_id'],
                how='left',
                suffixes=('', '_eff')
            )
        else:
            # 回退到只用genome
            effector_cols = ['genome', 'has_effector', 'effector_count',
                           'effector_types', 'effector_details']
            effector_cols = [c for c in effector_cols if c in effector_df.columns]
            result = result.merge(
                effector_df[effector_cols].drop_duplicates(subset=['genome']),
                on=['genome'],
                how='left',
                suffixes=('', '_eff')
            )

    logger.info(f"合并后记录数: {len(result)}")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Step 9: 综合报告生成",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-o', '--output', default='09_report', help='输出目录')
    parser.add_argument('--search', help='搜索结果文件')
    parser.add_argument('--motif', help='Motif过滤结果')
    parser.add_argument('--classified', help='分类结果')
    parser.add_argument('--rna', help='RNA筛选结果')
    parser.add_argument('--neighbors', help='邻近基因结果')
    parser.add_argument('--defense', help='防御岛分析结果')
    parser.add_argument('--effector', help='效应蛋白检测结果')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 9: 综合报告生成")
    logger.info("=" * 60)

    # 收集结果文件
    files = {
        'search': args.search,
        'motif': args.motif,
        'classified': args.classified,
        'rna': args.rna,
        'neighbors': args.neighbors,
        'defense': args.defense,
        'effector': args.effector
    }

    logger.info("\n加载结果文件...")
    result = merge_results(files, logger)

    if result is None or result.empty:
        logger.error("无法生成报告：没有结果数据")
        sys.exit(1)

    # 计算置信度评分
    logger.info("\n计算置信度评分...")
    scores = []
    evidences = []

    for _, row in result.iterrows():
        score, evidence = calculate_confidence_score(row)
        scores.append(score)
        evidences.append(evidence)

    result['confidence_score'] = scores
    result['confidence_evidence'] = evidences
    result['confidence_level'] = result['confidence_score'].apply(classify_confidence)

    # 按置信度排序
    result = result.sort_values('confidence_score', ascending=False)

    # 统计
    level_counts = result['confidence_level'].value_counts()
    logger.info("\n置信度分布:")
    for level in ['HIGH', 'MEDIUM', 'LOW', 'VERY_LOW']:
        count = level_counts.get(level, 0)
        pct = count / len(result) * 100
        logger.info(f"  {level}: {count} ({pct:.1f}%)")

    # 保存完整结果
    full_file = output_dir / "final_results.tsv"
    result.to_csv(full_file, sep='\t', index=False)
    logger.info(f"\n✓ 完整结果: {full_file}")

    # 保存高置信度候选
    high_conf = result[result['confidence_level'] == 'HIGH']
    if not high_conf.empty:
        high_file = output_dir / "high_confidence_candidates.tsv"
        high_conf.to_csv(high_file, sep='\t', index=False)
        logger.info(f"✓ 高置信度: {high_file} ({len(high_conf)} 条)")

    # 保存中等置信度候选
    medium_conf = result[result['confidence_level'] == 'MEDIUM']
    if not medium_conf.empty:
        medium_file = output_dir / "medium_confidence_candidates.tsv"
        medium_conf.to_csv(medium_file, sep='\t', index=False)
        logger.info(f"✓ 中等置信度: {medium_file} ({len(medium_conf)} 条)")

    # 生成摘要报告
    report_file = output_dir / "summary_report.md"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("# Retron系统挖掘报告\n\n")
        f.write(f"**生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("## 1. 流程概要\n\n")
        f.write("| 步骤 | 说明 |\n")
        f.write("|------|------|\n")
        f.write("| Step 1 | RT同源物搜索 |\n")
        f.write("| Step 2 | Motif过滤 (NAXXH/VTG) |\n")
        f.write("| Step 3 | RT分类与排除 |\n")
        f.write("| Step 4-5 | ncRNA结构预测 |\n")
        f.write("| Step 6-8 | 效应蛋白与防御岛分析 |\n\n")

        f.write("## 2. 结果统计\n\n")
        f.write(f"**总候选数**: {len(result)}\n\n")
        f.write("| 置信度 | 数量 | 比例 |\n")
        f.write("|--------|------|------|\n")
        for level in ['HIGH', 'MEDIUM', 'LOW', 'VERY_LOW']:
            count = level_counts.get(level, 0)
            pct = count / len(result) * 100
            f.write(f"| {level} | {count} | {pct:.1f}% |\n")

        f.write("\n## 3. 评分标准\n\n")
        f.write("| 维度 | 分值 | 条件 |\n")
        f.write("|------|------|------|\n")
        f.write("| Motif | 30 | NAXXH + VTG |\n")
        f.write("| Motif | 20 | NAXXH only |\n")
        f.write("| RNA茎长 | 15 | ≥20bp |\n")
        f.write("| RNA茎长 | 10 | ≥15bp |\n")
        f.write("| 分支G | 5 | 存在 |\n")
        f.write("| 效应蛋白 | 25 | 邻近有效应蛋白 |\n")
        f.write("| 防御岛 | 20 | 位于防御系统附近 |\n")

        f.write("\n## 4. 高置信度候选 (Top 20)\n\n")
        if not high_conf.empty:
            top20 = high_conf.head(20)
            f.write("| 基因组 | 评分 | 证据 |\n")
            f.write("|--------|------|------|\n")
            for _, row in top20.iterrows():
                f.write(f"| {row['genome']} | {row['confidence_score']} | {row['confidence_evidence']} |\n")
        else:
            f.write("*无高置信度候选*\n")

        f.write("\n## 5. 后续建议\n\n")
        f.write("1. **高置信度候选**: 可直接进行实验验证\n")
        f.write("2. **中等置信度**: 建议补充HHpred结构域分析\n")
        f.write("3. **低置信度**: 需要更多计算证据支持\n")

    logger.info(f"✓ 摘要报告: {report_file}")

    logger.info("\n" + "=" * 60)
    logger.info("报告生成完成!")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
