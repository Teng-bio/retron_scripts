#!/usr/bin/env python3
"""
Step 9: 综合报告生成与Retron类型分类 (Mestre et al., 2020 方法)

基于统计关联分析结果，对Retron候选进行分类和报告生成。

核心逻辑:
1. 整合各步骤结果 (RT候选、RNA结构、关联蛋白簇)
2. 基于关联蛋白的功能注释进行Retron类型预测
3. 生成最终候选列表和分类报告

与Millman方法的区别:
- Millman: 基于防御岛和效应蛋白的置信度评分
- Mestre: 基于统计关联蛋白的类型分类 (Type I-XIII)

输出:
  - final_candidates.tsv: 最终Retron候选列表
  - classification_report.tsv: 类型分类详情
  - summary_report.md: Markdown格式总结报告
"""

import argparse
import pandas as pd
import sys
import logging
from pathlib import Path
from datetime import datetime
from collections import defaultdict


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


# ============================================================
# Retron类型定义 (基于Mestre et al., 2020 Table 1)
# ============================================================
RETRON_TYPE_MARKERS = {
    'Type_I': {
        'markers': ['ATPase', 'AAA', 'P-loop', 'Toprim', 'PtuA'],
        'description': 'ATPase/Toprim伴随蛋白',
        'effector_type': 'ATPase',
    },
    'Type_II': {
        'markers': ['NDT', 'Nucleoside deoxyribosyltransferase', 'deoxyribosyltransferase'],
        'description': 'NDT核苷脱氧核糖基转移酶',
        'effector_type': 'NDT',
    },
    'Type_III': {
        'markers': ['PRTase', 'Phosphoribosyltransferase', 'PRPP'],
        'description': 'PRTase磷酸核糖基转移酶',
        'effector_type': 'PRTase',
    },
    'Type_IV': {
        'markers': ['HNH', 'endonuclease', 'homing'],
        'description': 'HNH核酸内切酶',
        'effector_type': 'HNH',
    },
    'Type_V': {
        'markers': ['DUF4297', 'DUF'],
        'description': 'DUF4297结构域蛋白',
        'effector_type': 'DUF4297',
    },
    'Type_VI': {
        'markers': ['transmembrane', '2TM', 'membrane protein', 'TM helix'],
        'description': '双跨膜蛋白',
        'effector_type': '2TM',
    },
    'Type_VII': {
        'markers': ['RNase H', 'DUF3800', 'RNaseH'],
        'description': 'RNase H2-like蛋白',
        'effector_type': 'RNaseH',
    },
    'Type_VIII': {
        'markers': ['OLD', 'RecB', 'nuclease', 'exonuclease'],
        'description': 'OLD/RecB核酸酶',
        'effector_type': 'OLD',
    },
    'Type_IX': {
        'markers': ['HEPN', 'RNase', 'abortive infection'],
        'description': 'HEPN结构域RNase',
        'effector_type': 'HEPN',
    },
    'Type_X': {
        'markers': ['cold shock', 'CspA', 'CSD'],
        'description': '冷休克蛋白',
        'effector_type': 'ColdShock',
    },
    'Type_XI': {
        'markers': ['ribosyltransferase', 'ADP-ribosyl', 'NAD'],
        'description': 'ADP核糖基转移酶',
        'effector_type': 'Ribosyltransferase',
    },
    'Type_XII': {
        'markers': ['toxin', 'MazF', 'RelE', 'ParE', 'antitoxin'],
        'description': '毒素-抗毒素系统',
        'effector_type': 'TA',
    },
    'Type_XIII': {
        'markers': ['zinc finger', 'SWIM', 'DNA-binding', 'HTH'],
        'description': 'DNA结合蛋白/锌指结构域',
        'effector_type': 'DNA-binding',
    },
}


def classify_by_associated_protein(product_annotation):
    """
    根据关联蛋白的功能注释预测Retron类型

    返回: (type_name, confidence)
    """
    if not product_annotation or not isinstance(product_annotation, str):
        return 'Unknown', 0

    product_lower = product_annotation.lower()
    matches = []

    for type_name, type_info in RETRON_TYPE_MARKERS.items():
        for marker in type_info['markers']:
            if marker.lower() in product_lower:
                matches.append((type_name, type_info))
                break

    if len(matches) == 1:
        return matches[0][0], 100  # 单一匹配，高置信度
    elif len(matches) > 1:
        return matches[0][0], 60  # 多重匹配，中等置信度
    else:
        return 'Unknown', 0


def load_results(files, logger):
    """加载各步骤结果"""
    dfs = {}

    for name, path in files.items():
        if path and Path(path).exists():
            try:
                dfs[name] = pd.read_csv(path, sep='\t')
                logger.info(f"  加载 {name}: {len(dfs[name])} 条")
            except Exception as e:
                logger.warning(f"  加载 {name} 失败: {e}")

    return dfs


def merge_results(dfs, logger):
    """
    合并各步骤结果

    主表: RT信息 (rt_info.tsv)
    关联: 显著关联蛋白簇、蛋白簇统计
    """
    # RT信息为主表
    if 'rt_info' not in dfs:
        logger.error("缺少RT信息表 (rt_info.tsv)")
        return None

    result = dfs['rt_info'].copy()
    logger.info(f"\nRT主表: {len(result)} 个RT")

    # 合并分类结果 (如果有)
    if 'classified' in dfs:
        classified = dfs['classified']
        # 尝试按genome合并
        if 'genome' in classified.columns:
            merge_cols = ['genome']
            if 'protein_id' in result.columns and 'protein_id' in classified.columns:
                merge_cols.append('protein_id')

            result = result.merge(
                classified.drop(columns=['contig'], errors='ignore').drop_duplicates(subset=merge_cols),
                on=merge_cols,
                how='left',
                suffixes=('', '_classified')
            )

    # 合并motif信息
    if 'motif' in dfs:
        motif = dfs['motif']
        if 'genome' in motif.columns:
            motif_cols = ['genome', 'has_naxxh', 'has_vtg', 'motif_classification']
            motif_cols = [c for c in motif_cols if c in motif.columns]
            result = result.merge(
                motif[motif_cols].drop_duplicates(subset=['genome']),
                on='genome',
                how='left',
                suffixes=('', '_motif')
            )

    # 合并RNA筛选结果
    if 'rna' in dfs:
        rna = dfs['rna']
        if 'genome' in rna.columns:
            rna_cols = ['genome', 'mfe', 'longest_stem', 'branch_g_count']
            rna_cols = [c for c in rna_cols if c in rna.columns]

            # 按genome分组，选择MFE最低的
            if 'mfe' in rna.columns:
                rna['mfe'] = pd.to_numeric(rna['mfe'], errors='coerce')
                best_rna = rna.loc[rna.groupby('genome')['mfe'].idxmin()]
                result = result.merge(
                    best_rna[rna_cols].drop_duplicates(subset=['genome']),
                    on='genome',
                    how='left',
                    suffixes=('', '_rna')
                )

    return result


def assign_associated_clusters(result_df, cluster_stats_df, matrix_df, neighbor_matrix_df, logger):
    """
    为每个RT分配显著关联的蛋白簇

    逻辑:
    1. 从共现矩阵找出每个RT有哪些蛋白簇
    2. 筛选那些是显著关联的簇
    3. 获取簇的功能注释 (从neighbor_matrix)
    """
    if matrix_df is None or cluster_stats_df is None:
        logger.warning("缺少聚类数据，跳过关联蛋白分配")
        return result_df

    # 获取显著关联的簇列表
    significant_clusters = set()
    if 'is_significant' in cluster_stats_df.columns:
        significant_clusters = set(
            cluster_stats_df[cluster_stats_df['is_significant'] == True]['cluster'].tolist()
        )
    elif 'phyvalue' in cluster_stats_df.columns:
        significant_clusters = set(
            cluster_stats_df[cluster_stats_df['phyvalue'] >= 2]['cluster'].tolist()
        )

    logger.info(f"显著关联蛋白簇数: {len(significant_clusters)}")

    # 构建簇 -> 功能注释的映射 (从neighbor_matrix获取代表性注释)
    cluster_annotations = {}
    if neighbor_matrix_df is not None and 'product' in neighbor_matrix_df.columns:
        for cluster in significant_clusters:
            # 找到属于该簇的邻近基因
            cluster_members = neighbor_matrix_df[
                neighbor_matrix_df['neighbor_id'].str.contains(f"|{cluster}|", regex=False, na=False)
            ]
            if not cluster_members.empty:
                # 取最常见的非空product注释
                products = cluster_members['product'].dropna().tolist()
                if products:
                    from collections import Counter
                    product_counts = Counter(products)
                    most_common = product_counts.most_common(1)[0][0]
                    cluster_annotations[cluster] = most_common

    # 为每个RT分配关联蛋白簇
    associated_clusters_list = []
    predicted_types_list = []
    type_confidence_list = []

    for rt_id in result_df['rt_id']:
        if rt_id not in matrix_df.index:
            associated_clusters_list.append('')
            predicted_types_list.append('Unknown')
            type_confidence_list.append(0)
            continue

        # 找出该RT有的显著关联簇
        rt_clusters = []
        for cluster in significant_clusters:
            if cluster in matrix_df.columns and matrix_df.loc[rt_id, cluster] == 1:
                rt_clusters.append(cluster)

        # 获取注释并预测类型
        annotations = [cluster_annotations.get(c, '') for c in rt_clusters if c in cluster_annotations]
        all_annotations = ' '.join(annotations)

        # 预测Retron类型
        predicted_type, confidence = classify_by_associated_protein(all_annotations)

        associated_clusters_list.append(';'.join(rt_clusters[:5]))  # 最多显示5个
        predicted_types_list.append(predicted_type)
        type_confidence_list.append(confidence)

    result_df['associated_clusters'] = associated_clusters_list
    result_df['predicted_type'] = predicted_types_list
    result_df['type_confidence'] = type_confidence_list

    # 统计类型分布
    type_counts = result_df['predicted_type'].value_counts()
    logger.info(f"\n预测类型分布:")
    for ptype, count in type_counts.items():
        logger.info(f"  {ptype}: {count}")

    return result_df


def calculate_final_score(row):
    """
    计算最终综合评分 (基于Mestre方法)

    评分维度:
    - Motif (30分): NAXXH + VTG
    - RNA结构 (25分): 茎长度、分支G
    - 关联蛋白 (30分): 有显著关联蛋白簇
    - 类型预测 (15分): 能预测到具体类型
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
    elif row.get('has_naxxh', False):
        score += 15
        evidence.append("Motif:NAXXH")
    elif row.get('has_vtg', False):
        score += 10
        evidence.append("Motif:VTG")

    # 2. RNA结构评分 (25分)
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

    branch_g = row.get('branch_g_count', 0)
    if pd.notna(branch_g) and int(branch_g) > 0:
        score += 5
        evidence.append(f"BranchG:{int(branch_g)}")

    mfe = row.get('mfe', 0)
    if pd.notna(mfe):
        mfe = float(mfe)
        if mfe < -25:
            score += 5
        elif mfe < -20:
            score += 3

    # 3. 关联蛋白评分 (30分)
    assoc_clusters = row.get('associated_clusters', '')
    if assoc_clusters and assoc_clusters != '':
        n_clusters = len(assoc_clusters.split(';'))
        if n_clusters >= 3:
            score += 30
            evidence.append(f"AssocProt:{n_clusters}")
        elif n_clusters >= 1:
            score += 20
            evidence.append(f"AssocProt:{n_clusters}")

    # 4. 类型预测评分 (15分)
    predicted_type = row.get('predicted_type', 'Unknown')
    type_conf = row.get('type_confidence', 0)
    if predicted_type != 'Unknown':
        if type_conf >= 80:
            score += 15
            evidence.append(f"Type:{predicted_type}")
        elif type_conf >= 50:
            score += 10
            evidence.append(f"Type:{predicted_type}?")

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


def generate_summary_report(result_df, output_dir, logger):
    """生成Markdown格式总结报告"""
    report_file = output_dir / "summary_report.md"

    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("# Retron系统挖掘报告 (Mestre方法)\n\n")
        f.write(f"**生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("## 1. 方法说明\n\n")
        f.write("本报告采用 Mestre et al., 2020 (NAR) 的统计关联方法:\n")
        f.write("- 提取RT上下游±30kb内的所有ORF\n")
        f.write("- 使用MMseqs2进行蛋白序列聚类\n")
        f.write("- 通过Phyvalue统计检验筛选显著关联蛋白簇\n")
        f.write("- 基于关联蛋白功能预测Retron类型 (Type I-XIII)\n\n")

        f.write("## 2. 流程概要\n\n")
        f.write("| 步骤 | 说明 |\n")
        f.write("|------|------|\n")
        f.write("| Step 1-3 | RT搜索与Motif筛选 |\n")
        f.write("| Step 4-5 | ncRNA结构预测 |\n")
        f.write("| Step 6 | 邻近蛋白提取 (±30kb) |\n")
        f.write("| Step 7 | MMseqs2聚类 |\n")
        f.write("| Step 8 | Phyvalue统计分析 |\n")
        f.write("| Step 9 | 类型分类与报告 |\n\n")

        f.write("## 3. 结果统计\n\n")
        f.write(f"**总候选数**: {len(result_df)}\n\n")

        # 置信度分布
        if 'confidence_level' in result_df.columns:
            level_counts = result_df['confidence_level'].value_counts()
            f.write("### 3.1 置信度分布\n\n")
            f.write("| 置信度 | 数量 | 比例 |\n")
            f.write("|--------|------|------|\n")
            for level in ['HIGH', 'MEDIUM', 'LOW', 'VERY_LOW']:
                count = level_counts.get(level, 0)
                pct = count / len(result_df) * 100
                f.write(f"| {level} | {count} | {pct:.1f}% |\n")
            f.write("\n")

        # 类型分布
        if 'predicted_type' in result_df.columns:
            type_counts = result_df['predicted_type'].value_counts()
            f.write("### 3.2 Retron类型分布\n\n")
            f.write("| 类型 | 数量 | 描述 |\n")
            f.write("|------|------|------|\n")
            for ptype, count in type_counts.items():
                desc = RETRON_TYPE_MARKERS.get(ptype, {}).get('description', '-')
                f.write(f"| {ptype} | {count} | {desc} |\n")
            f.write("\n")

        f.write("## 4. 评分标准\n\n")
        f.write("| 维度 | 分值 | 条件 |\n")
        f.write("|------|------|------|\n")
        f.write("| Motif | 30 | NAXXH + VTG |\n")
        f.write("| Motif | 20 | NAXXH only |\n")
        f.write("| RNA茎长 | 15 | ≥20bp |\n")
        f.write("| RNA茎长 | 10 | ≥15bp |\n")
        f.write("| 分支G | 5 | 存在 |\n")
        f.write("| 关联蛋白 | 30 | ≥3个显著关联簇 |\n")
        f.write("| 关联蛋白 | 20 | 1-2个显著关联簇 |\n")
        f.write("| 类型预测 | 15 | 高置信度类型预测 |\n\n")

        # Top候选
        if 'confidence_score' in result_df.columns:
            high_conf = result_df[result_df['confidence_level'] == 'HIGH']
            f.write("## 5. 高置信度候选 (Top 20)\n\n")
            if not high_conf.empty:
                f.write("| RT ID | 基因组 | 类型 | 评分 | 证据 |\n")
                f.write("|-------|--------|------|------|------|\n")
                for _, row in high_conf.head(20).iterrows():
                    rt_id = row.get('rt_id', '-')[:30]
                    genome = row.get('genome', '-')[:20]
                    ptype = row.get('predicted_type', '-')
                    score = row.get('confidence_score', 0)
                    evidence = row.get('confidence_evidence', '-')[:40]
                    f.write(f"| {rt_id} | {genome} | {ptype} | {score} | {evidence} |\n")
            else:
                f.write("*无高置信度候选*\n")
            f.write("\n")

        f.write("## 6. 后续建议\n\n")
        f.write("1. **高置信度候选**: 优先进行实验验证\n")
        f.write("2. **有类型预测的候选**: 参考对应类型的已知Retron进行功能验证\n")
        f.write("3. **Unknown类型**: 可能是新型Retron，值得深入研究\n")
        f.write("4. **关联蛋白注释**: 对显著关联蛋白簇进行InterProScan/HHpred注释\n")

    logger.info(f"✓ 总结报告: {report_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Step 9: 综合报告生成与Retron类型分类 (Mestre方法)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-o', '--output', default='09_report', help='输出目录')

    # 输入文件
    parser.add_argument('--rt-info', help='RT信息表 (step06/rt_info.tsv)')
    parser.add_argument('--motif', help='Motif过滤结果')
    parser.add_argument('--classified', help='RT分类结果')
    parser.add_argument('--rna', help='RNA筛选结果')
    parser.add_argument('--cluster-stats', help='Phyvalue分析结果 (step08/phyvalue_analysis.tsv)')
    parser.add_argument('--matrix', help='RT-簇共现矩阵 (step07/rt_cluster_matrix.tsv)')
    parser.add_argument('--neighbor-matrix', help='邻近基因矩阵 (step06/neighbor_matrix.tsv)')

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 9: 综合报告生成与Retron类型分类 (Mestre方法)")
    logger.info("=" * 60)

    # 收集结果文件
    files = {
        'rt_info': args.rt_info,
        'motif': args.motif,
        'classified': args.classified,
        'rna': args.rna,
    }

    logger.info("\n加载结果文件...")
    dfs = load_results(files, logger)

    if 'rt_info' not in dfs:
        logger.error("缺少RT信息表，无法生成报告")
        sys.exit(1)

    # 合并结果
    result = merge_results(dfs, logger)

    if result is None or result.empty:
        logger.error("无法生成报告：没有结果数据")
        sys.exit(1)

    # 加载聚类和关联分析结果
    matrix_df = None
    cluster_stats_df = None
    neighbor_matrix_df = None

    if args.matrix and Path(args.matrix).exists():
        matrix_df = pd.read_csv(args.matrix, sep='\t', index_col=0)
        logger.info(f"  加载共现矩阵: {matrix_df.shape}")

    if args.cluster_stats and Path(args.cluster_stats).exists():
        cluster_stats_df = pd.read_csv(args.cluster_stats, sep='\t')
        logger.info(f"  加载Phyvalue分析: {len(cluster_stats_df)} 条")

    if args.neighbor_matrix and Path(args.neighbor_matrix).exists():
        neighbor_matrix_df = pd.read_csv(args.neighbor_matrix, sep='\t')
        logger.info(f"  加载邻近基因矩阵: {len(neighbor_matrix_df)} 条")

    # 分配关联蛋白簇并预测类型
    if matrix_df is not None and cluster_stats_df is not None:
        logger.info("\n分配关联蛋白簇...")
        result = assign_associated_clusters(
            result, cluster_stats_df, matrix_df, neighbor_matrix_df, logger
        )

    # 计算综合评分
    logger.info("\n计算综合评分...")
    scores = []
    evidences = []

    for _, row in result.iterrows():
        score, evidence = calculate_final_score(row)
        scores.append(score)
        evidences.append(evidence)

    result['confidence_score'] = scores
    result['confidence_evidence'] = evidences
    result['confidence_level'] = result['confidence_score'].apply(classify_confidence)

    # 按评分排序
    result = result.sort_values('confidence_score', ascending=False)

    # 统计
    level_counts = result['confidence_level'].value_counts()
    logger.info("\n置信度分布:")
    for level in ['HIGH', 'MEDIUM', 'LOW', 'VERY_LOW']:
        count = level_counts.get(level, 0)
        pct = count / len(result) * 100
        logger.info(f"  {level}: {count} ({pct:.1f}%)")

    # 保存结果
    full_file = output_dir / "final_candidates.tsv"
    result.to_csv(full_file, sep='\t', index=False)
    logger.info(f"\n✓ 完整结果: {full_file}")

    # 高置信度候选
    high_conf = result[result['confidence_level'] == 'HIGH']
    if not high_conf.empty:
        high_file = output_dir / "high_confidence_candidates.tsv"
        high_conf.to_csv(high_file, sep='\t', index=False)
        logger.info(f"✓ 高置信度: {high_file} ({len(high_conf)} 条)")

    # 按类型分组
    if 'predicted_type' in result.columns:
        for ptype in result['predicted_type'].unique():
            if ptype == 'Unknown':
                continue
            type_df = result[result['predicted_type'] == ptype]
            if not type_df.empty:
                type_file = output_dir / f"candidates_{ptype}.tsv"
                type_df.to_csv(type_file, sep='\t', index=False)
                logger.info(f"✓ {ptype}: {type_file} ({len(type_df)} 条)")

    # 生成总结报告
    generate_summary_report(result, output_dir, logger)

    # 分类统计
    class_file = output_dir / "classification_summary.tsv"
    if 'predicted_type' in result.columns:
        type_summary = result.groupby('predicted_type').agg({
            'rt_id': 'count',
            'confidence_score': ['mean', 'max'],
            'genome': 'nunique'
        }).round(1)
        type_summary.columns = ['count', 'avg_score', 'max_score', 'unique_genomes']
        type_summary = type_summary.sort_values('count', ascending=False)
        type_summary.to_csv(class_file, sep='\t')
        logger.info(f"✓ 分类汇总: {class_file}")

    logger.info("\n" + "=" * 60)
    logger.info("报告生成完成!")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
