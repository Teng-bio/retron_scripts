#!/usr/bin/env python3
"""
Step 2: Retron RT Motif过滤

检测RT序列中的Retron特征性motif：
- Region X: NAXXH motif (区分Retron RT与其他RT的关键)
- Region Y: VTG motif (C末端特征)

参考: Millman et al., 2020 Cell
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


def filter_by_motifs(df, sequences, logger):
    """
    基于motif过滤RT候选

    返回: 添加了motif信息的DataFrame
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

        # 检测motifs
        naxxh_hits = find_naxxh_motif(sequence)
        vtg_hits = find_vtg_motif(sequence)

        has_naxxh = len(naxxh_hits) > 0
        has_vtg = len(vtg_hits) > 0

        # 分类
        if has_naxxh and has_vtg:
            classification = 'high_confidence'
        elif has_naxxh:
            classification = 'medium_confidence'
        elif has_vtg:
            classification = 'low_confidence'
        else:
            classification = 'no_motif'

        result = row.to_dict()
        result.update({
            'rt_length': len(sequence),
            'has_naxxh': has_naxxh,
            'naxxh_motifs': ';'.join([f"{m[2]}@{m[0]}" for m in naxxh_hits]) if naxxh_hits else '',
            'has_vtg': has_vtg,
            'vtg_motifs': ';'.join([f"{m[2]}@{m[0]}" for m in vtg_hits]) if vtg_hits else '',
            'motif_classification': classification,
            'length_ok': length_ok
        })
        results.append(result)

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Retron RT Motif过滤 (NAXXH + VTG)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Step1的搜索结果 (complete_results.tsv)')
    parser.add_argument('-f', '--faa-dir', required=True,
                        help='蛋白序列目录 (.faa文件)')
    parser.add_argument('-o', '--output', default='02_motif_filtered',
                        help='输出目录')
    parser.add_argument('--keep-all', action='store_true',
                        help='保留所有结果（包括无motif的）')
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
    logger.info("Step 2: Retron RT Motif过滤")
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

    # Motif过滤
    logger.info("\n检测Retron特征motif...")
    result_df = filter_by_motifs(df, sequences, logger)

    # 统计
    stats = result_df['motif_classification'].value_counts()
    logger.info("\nMotif分类统计:")
    for cls, count in stats.items():
        logger.info(f"  {cls}: {count}")

    # 保存结果
    # 完整结果
    full_output = output_dir / "motif_analysis_full.tsv"
    result_df.to_csv(full_output, sep='\t', index=False)
    logger.info(f"\n✓ 完整结果: {full_output}")

    # 过滤后结果（至少有一个motif）
    if not args.keep_all:
        filtered_df = result_df[result_df['motif_classification'] != 'no_motif']
    else:
        filtered_df = result_df

    filtered_output = output_dir / "motif_filtered.tsv"
    filtered_df.to_csv(filtered_output, sep='\t', index=False)
    logger.info(f"✓ 过滤结果: {filtered_output} ({len(filtered_df)} 条)")

    # 高置信度结果
    high_conf_df = result_df[result_df['motif_classification'] == 'high_confidence']
    if not high_conf_df.empty:
        high_conf_output = output_dir / "high_confidence.tsv"
        high_conf_df.to_csv(high_conf_output, sep='\t', index=False)
        logger.info(f"✓ 高置信度: {high_conf_output} ({len(high_conf_df)} 条)")

    # 统计报告
    stats_file = output_dir / "motif_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("Retron RT Motif过滤统计\n")
        f.write("=" * 40 + "\n")
        f.write(f"输入: {len(df)} 条\n")
        f.write(f"检测到序列: {len(result_df)} 条\n\n")
        f.write("分类统计:\n")
        for cls, count in stats.items():
            pct = count / len(result_df) * 100
            f.write(f"  {cls}: {count} ({pct:.1f}%)\n")
        f.write(f"\n过滤后保留: {len(filtered_df)} 条\n")

    logger.info("\n" + "=" * 60)
    logger.info(f"过滤完成! {len(df)} -> {len(filtered_df)} 条")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
