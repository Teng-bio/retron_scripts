#!/usr/bin/env python3
"""
Step 2: Retron RT Motif Tier 分层分类

基于保守 motif 的软分类系统（不进行硬过滤）：
- Tier 1 (高置信度): NAXXH + VTG 都存在 - 经典 Retron RT 标志
- Tier 2 (潜在新型): 部分 motif (只有NAXXH或VTG) - 可能是新型 Retron
- Tier 3 (排除): 有 YADD 标记 (Group II Intron) 或无任何 Retron motif

检测的 motif：
- NAXXH (Region X): [NQDE]A..[HY] - 区分 Retron RT 与其他 RT 的关键
- VTG (Region Y): V[TS]G - C末端特征
- YADD: Y[AV]DD - Group II Intron RT 标志 (排除标记)

参考:
- Millman et al., 2020 Cell (Retron motifs)
- Mestre et al., 2020 NAR (统计关联方法)
"""

import argparse
import pandas as pd
import re
import sys
import logging
from pathlib import Path
from Bio import SeqIO


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


def find_naxxh_motif(sequence):
    """
    检测NAXXH motif (Region X)

    NAXXH: N=任意, A=丙氨酸, X=任意, X=任意, H=组氨酸
    允许一定变异: [NQDE]A..[HY]

    返回: [(start, end, motif_seq), ...]
    """
    # 严格模式: NA..H
    strict_pattern = r'NA..[HY]'
    # 宽松模式: [NQDE]A..[HY]
    relaxed_pattern = r'[NQDE]A..[HY]'

    matches = []

    # 先找严格匹配
    for m in re.finditer(strict_pattern, sequence):
        matches.append((m.start(), m.end(), m.group(), 'strict'))

    # 如果没有严格匹配，找宽松匹配
    if not matches:
        for m in re.finditer(relaxed_pattern, sequence):
            matches.append((m.start(), m.end(), m.group(), 'relaxed'))

    return matches


def find_vtg_motif(sequence):
    """
    检测VTG motif (Region Y, C末端)

    VTG: 缬氨酸-苏氨酸-甘氨酸
    通常位于C末端100aa内

    返回: [(start, end, motif_seq), ...]
    """
    # VTG及其变体
    pattern = r'V[TS]G'

    matches = []
    seq_len = len(sequence)

    # 只在C末端150aa范围内搜索
    search_region = sequence[-150:] if seq_len > 150 else sequence
    offset = max(0, seq_len - 150)

    for m in re.finditer(pattern, search_region):
        matches.append((m.start() + offset, m.end() + offset, m.group()))

    return matches


def find_yadd_motif(sequence):
    """
    检测YADD motif - Group II Intron RT 的特征标记

    YADD: Y[AV]DD - Group II Intron RT 活性位点
    如果存在此 motif，很可能是 Group II Intron RT 而非 Retron RT

    返回: [(start, end, motif_seq), ...]
    """
    # YADD 及其变体
    pattern = r'Y[AV]DD'

    matches = []
    for m in re.finditer(pattern, sequence):
        matches.append((m.start(), m.end(), m.group()))

    return matches


def check_rt_length(sequence, min_len=250, max_len=600):
    """检查RT长度是否在合理范围内"""
    return min_len <= len(sequence) <= max_len


def extract_core_number(name):
    """提取核心数字用于文件匹配"""
    match = re.search(r'(\d{9,}(?:\.\d+)?)', str(name))
    return match.group(1) if match else str(name)


def load_sequences_from_faa(faa_dir, genome_list, logger):
    """从faa文件加载指定基因组的RT序列（与step01兼容）"""
    faa_dir = Path(faa_dir)
    sequences = {}

    # 构建文件索引（支持多种匹配方式）
    # 使用rglob递归搜索子目录中的faa文件
    faa_files = {}

    # 先尝试顶层目录
    for f in faa_dir.glob('*.faa'):
        faa_files[f.stem] = f
        core = extract_core_number(f.stem)
        if core not in faa_files:
            faa_files[core] = f

    # 如果顶层没有找到，递归搜索子目录
    if not faa_files:
        logger.info("顶层目录未找到faa文件，递归搜索子目录...")
        for f in faa_dir.rglob('*.faa'):
            faa_files[f.stem] = f
            # 也用父目录名作为索引（如 003C31/003C31.faa -> 003C31）
            faa_files[f.parent.name] = f
            core = extract_core_number(f.stem)
            if core not in faa_files:
                faa_files[core] = f

    logger.info(f"找到 {len(faa_files)} 个faa文件索引")

    loaded_genomes = set()
    for genome in genome_list:
        # 查找对应的faa文件
        faa_file = faa_files.get(genome)
        if not faa_file:
            core = extract_core_number(genome)
            faa_file = faa_files.get(core)

        if faa_file and faa_file.exists():
            genome_name = faa_file.stem  # 使用文件名作为genome_name
            if genome_name in loaded_genomes:
                continue
            loaded_genomes.add(genome_name)

            for record in SeqIO.parse(faa_file, 'fasta'):
                # 使用与step01相同的格式: genome_name___seq_id
                key = f"{genome_name}___{record.id}"
                sequences[key] = str(record.seq)
                # 也存储原始ID以便直接匹配
                sequences[record.id] = str(record.seq)

    logger.info(f"加载了 {len(sequences)} 条蛋白序列 (来自 {len(loaded_genomes)} 个基因组)")
    return sequences


def classify_by_motifs(df, sequences, logger):
    """
    基于 motif 的 Tier 分层分类

    Tier 1: NAXXH + VTG (且无 YADD) - 高置信度 Retron RT
    Tier 2: 只有 NAXXH 或只有 VTG (且无 YADD) - 潜在新型 Retron
    Tier 3: 有 YADD 标记 或 无任何 Retron motif - 排除候选

    返回: 添加了 tier 信息的 DataFrame
    """
    results = []
    found_count = 0
    not_found_count = 0

    for _, row in df.iterrows():
        genome = row['genome']
        contig = row.get('contig', '')

        # 尝试多种方式获取序列（与step01的genome___seq_id格式兼容）
        sequence = None
        tried_keys = []

        # 方式1: genome___contig
        key1 = f"{genome}___{contig}"
        tried_keys.append(key1)
        sequence = sequences.get(key1)

        # 方式2: 直接用contig（蛋白质ID）
        if not sequence:
            tried_keys.append(contig)
            sequence = sequences.get(contig)

        # 方式3: 如果contig包含genome___前缀，去掉后再试
        if not sequence and '___' in contig:
            parts = contig.split('___', 1)
            if len(parts) == 2:
                tried_keys.append(parts[1])
                sequence = sequences.get(parts[1])

        if not sequence:
            not_found_count += 1
            if not_found_count <= 5:
                logger.debug(f"未找到序列: 尝试了 {tried_keys}")
            continue

        found_count += 1

        # 检查长度
        length_ok = check_rt_length(sequence)

        # 检测所有 motifs
        naxxh_hits = find_naxxh_motif(sequence)
        vtg_hits = find_vtg_motif(sequence)
        yadd_hits = find_yadd_motif(sequence)

        has_naxxh = len(naxxh_hits) > 0
        has_vtg = len(vtg_hits) > 0
        has_yadd = len(yadd_hits) > 0

        # Tier 分层分类
        if has_yadd:
            # 有 Group II Intron 标记 -> Tier 3 (排除)
            tier = 3
            tier_reason = 'YADD_marker'
            recommended_action = 'exclude'
        elif has_naxxh and has_vtg:
            # 经典 Retron RT 标志 -> Tier 1 (高置信度)
            tier = 1
            tier_reason = 'NAXXH+VTG'
            recommended_action = 'keep'
        elif has_naxxh:
            # 只有 NAXXH -> Tier 2 (潜在新型)
            tier = 2
            tier_reason = 'NAXXH_only'
            recommended_action = 'review'
        elif has_vtg:
            # 只有 VTG -> Tier 2 (潜在新型)
            tier = 2
            tier_reason = 'VTG_only'
            recommended_action = 'review'
        else:
            # 无任何 Retron motif -> Tier 3 (排除)
            tier = 3
            tier_reason = 'no_retron_motif'
            recommended_action = 'exclude'

        result = row.to_dict()
        result.update({
            'rt_length': len(sequence),
            # Motif 检测结果
            'has_naxxh': has_naxxh,
            'naxxh_motifs': ';'.join([f"{m[2]}@{m[0]}" for m in naxxh_hits]) if naxxh_hits else '',
            'has_vtg': has_vtg,
            'vtg_motifs': ';'.join([f"{m[2]}@{m[0]}" for m in vtg_hits]) if vtg_hits else '',
            'has_yadd': has_yadd,
            'yadd_motifs': ';'.join([f"{m[2]}@{m[0]}" for m in yadd_hits]) if yadd_hits else '',
            # Tier 分类
            'tier': tier,
            'tier_reason': tier_reason,
            'recommended_action': recommended_action,
            'length_ok': length_ok
        })
        results.append(result)

    logger.info(f"成功获取序列: {found_count}, 未找到: {not_found_count}")
    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Retron RT Motif Tier 分层分类",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Tier 分类说明:
  Tier 1: NAXXH + VTG (高置信度) - 经典 Retron RT 标志
  Tier 2: 只有 NAXXH 或 VTG (潜在新型) - 需人工审查
  Tier 3: 有 YADD 标记或无 motif (排除) - 可能是 Group II Intron 或其他 RT

输出文件:
  motif_tier_all.tsv    - 所有候选的完整 tier 分类结果
  tier1_candidates.tsv  - Tier 1 高置信度候选 (推荐保留)
  tier2_candidates.tsv  - Tier 2 潜在新型候选 (需审查)
  tier3_excluded.tsv    - Tier 3 排除候选 (建议排除)
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Step1的搜索结果 (complete_results.tsv)')
    parser.add_argument('-f', '--faa-dir', required=True,
                        help='蛋白序列目录 (.faa文件)')
    parser.add_argument('-o', '--output', default='02_motif_filtered',
                        help='输出目录')
    parser.add_argument('--exclude-tier3', action='store_true',
                        help='输出时排除 Tier 3 候选 (默认保留所有)')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    # 验证输入
    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"输入文件不存在: {input_file}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 2: Retron RT Motif Tier 分层分类")
    logger.info("=" * 60)

    # 读取搜索结果
    df = pd.read_csv(input_file, sep='\t')
    logger.info(f"读取到 {len(df)} 条RT候选")

    # 获取唯一基因组列表
    genomes = df['genome'].unique().tolist()
    logger.info(f"涉及 {len(genomes)} 个基因组")

    # 加载序列
    logger.info("\n加载蛋白序列...")
    sequences = load_sequences_from_faa(args.faa_dir, genomes, logger)

    # Tier 分层分类
    logger.info("\n检测 Retron 特征 motif 并进行 Tier 分类...")
    result_df = classify_by_motifs(df, sequences, logger)

    if result_df.empty:
        logger.warning("未能获取任何序列进行分类")
        sys.exit(1)

    # Tier 统计
    tier_stats = result_df['tier'].value_counts().sort_index()
    reason_stats = result_df['tier_reason'].value_counts()

    logger.info("\n" + "=" * 40)
    logger.info("Tier 分类统计:")
    logger.info("=" * 40)
    for tier_num in [1, 2, 3]:
        count = tier_stats.get(tier_num, 0)
        pct = count / len(result_df) * 100 if len(result_df) > 0 else 0
        tier_names = {1: '高置信度', 2: '潜在新型', 3: '排除候选'}
        logger.info(f"  Tier {tier_num} ({tier_names[tier_num]}): {count} ({pct:.1f}%)")

    logger.info("\n分类原因统计:")
    for reason, count in reason_stats.items():
        logger.info(f"  {reason}: {count}")

    # 保存结果
    # 1. 完整结果（所有 tier）
    all_output = output_dir / "motif_tier_all.tsv"
    result_df.to_csv(all_output, sep='\t', index=False)
    logger.info(f"\n✓ 完整结果: {all_output}")

    # 2. 分 tier 输出
    tier1_df = result_df[result_df['tier'] == 1]
    tier2_df = result_df[result_df['tier'] == 2]
    tier3_df = result_df[result_df['tier'] == 3]

    if not tier1_df.empty:
        tier1_output = output_dir / "tier1_candidates.tsv"
        tier1_df.to_csv(tier1_output, sep='\t', index=False)
        logger.info(f"✓ Tier 1 (高置信度): {tier1_output} ({len(tier1_df)} 条)")

    if not tier2_df.empty:
        tier2_output = output_dir / "tier2_candidates.tsv"
        tier2_df.to_csv(tier2_output, sep='\t', index=False)
        logger.info(f"✓ Tier 2 (潜在新型): {tier2_output} ({len(tier2_df)} 条)")

    if not tier3_df.empty:
        tier3_output = output_dir / "tier3_excluded.tsv"
        tier3_df.to_csv(tier3_output, sep='\t', index=False)
        logger.info(f"✓ Tier 3 (排除): {tier3_output} ({len(tier3_df)} 条)")

    # 3. 兼容旧流程：输出 motif_filtered.tsv (Tier 1 + Tier 2)
    if args.exclude_tier3:
        filtered_df = result_df[result_df['tier'] <= 2]
    else:
        filtered_df = result_df

    filtered_output = output_dir / "motif_filtered.tsv"
    filtered_df.to_csv(filtered_output, sep='\t', index=False)
    logger.info(f"✓ 过滤结果: {filtered_output} ({len(filtered_df)} 条)")

    # 4. 高置信度输出（兼容旧流程）
    if not tier1_df.empty:
        high_conf_output = output_dir / "high_confidence.tsv"
        tier1_df.to_csv(high_conf_output, sep='\t', index=False)

    # 统计报告
    stats_file = output_dir / "tier_classification_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("Retron RT Motif Tier 分层分类报告\n")
        f.write("=" * 50 + "\n\n")

        f.write("Tier 分类标准:\n")
        f.write("  Tier 1: NAXXH + VTG (高置信度) - 经典 Retron RT\n")
        f.write("  Tier 2: 只有 NAXXH 或 VTG (潜在新型) - 需审查\n")
        f.write("  Tier 3: YADD 标记或无 motif (排除) - 可能非 Retron\n\n")

        f.write(f"输入: {len(df)} 条候选\n")
        f.write(f"成功分析: {len(result_df)} 条\n\n")

        f.write("Tier 分布:\n")
        for tier_num in [1, 2, 3]:
            count = tier_stats.get(tier_num, 0)
            pct = count / len(result_df) * 100 if len(result_df) > 0 else 0
            tier_names = {1: '高置信度', 2: '潜在新型', 3: '排除候选'}
            f.write(f"  Tier {tier_num} ({tier_names[tier_num]}): {count} ({pct:.1f}%)\n")

        f.write("\n详细分类原因:\n")
        for reason, count in reason_stats.items():
            f.write(f"  {reason}: {count}\n")

        f.write(f"\n推荐操作:\n")
        f.write(f"  - Tier 1: 直接进入下一步 ({len(tier1_df)} 条)\n")
        f.write(f"  - Tier 2: 人工审查后决定 ({len(tier2_df)} 条)\n")
        f.write(f"  - Tier 3: 建议排除 ({len(tier3_df)} 条)\n")

    logger.info(f"✓ 统计报告: {stats_file}")

    # 打印总结
    logger.info("\n" + "=" * 60)
    logger.info("分类完成!")
    logger.info("=" * 60)
    logger.info(f"Tier 1 (高置信度): {len(tier1_df)} 条 -> 推荐保留")
    logger.info(f"Tier 2 (潜在新型): {len(tier2_df)} 条 -> 需审查")
    logger.info(f"Tier 3 (排除):     {len(tier3_df)} 条 -> 建议排除")
    logger.info(f"\n下一步: 使用 motif_filtered.tsv 或 tier1_candidates.tsv 进入 Step 3")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
