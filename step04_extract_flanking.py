#!/usr/bin/env python3
"""
Step 2: 上游盲搜序列提取 (Upstream Blind Extraction)

针对黏细菌(High-GC)的Retron挖掘优化版本。
根据文献(Inouye et al., 1989; Rice et al., 1995)，Retron的ncRNA(msr-msd)
严格位于RT基因的5'端上游，且有时会延伸覆盖到RT基因的起始部分。

核心策略：
- 只提取上游序列，不提取下游（避免虚假折叠）
- 包含RT头部以覆盖潜在的msd-RT重叠区
- 不依赖基因注释，直接基于RT坐标盲切

提取范围：
- 正链(+): [RT_Start - upstream, RT_Start + head_overlap]
- 负链(-): [RT_End - head_overlap, RT_End + upstream]，然后反向互补

输入：Step1的complete_results.tsv + 基因组fasta文件
输出：upstream_sequences.fasta
"""

import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import logging
import sys


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


def extract_core_number(name):
    """提取核心数字用于文件匹配"""
    match = re.search(r'(\d{9,}(?:\.\d+)?)', str(name))
    return match.group(1) if match else str(name)


def build_fasta_index(fasta_dir, logger):
    """构建fasta文件索引"""
    fasta_dir = Path(fasta_dir)
    index = {}
    for ext in ['.fasta', '.fna', '.fa']:
        for f in fasta_dir.glob(f'*{ext}'):
            index[f.stem] = f
            core = extract_core_number(f.stem)
            if core not in index:
                index[core] = f
    logger.info(f"索引了 {len(set(index.values()))} 个fasta文件")
    return index


def find_fasta_file(genome_name, fasta_index):
    """根据genome_name查找对应的fasta文件"""
    if genome_name in fasta_index:
        return fasta_index[genome_name]
    core = extract_core_number(genome_name)
    return fasta_index.get(core)


def load_genome_sequences(fasta_file):
    """加载基因组序列到字典"""
    sequences = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences[record.id] = record
        if '.' in record.id:
            sequences[record.id.split('.')[0]] = record
    return sequences


def find_contig(genome_seqs, contig_name):
    """查找contig，支持多种名称变体"""
    search_names = [contig_name]
    if contig_name.startswith('NZ_'):
        search_names.append(contig_name[3:])
    else:
        search_names.append('NZ_' + contig_name)
    if '.' in contig_name:
        search_names.append(contig_name.split('.')[0])

    for name in search_names:
        if name in genome_seqs:
            return genome_seqs[name]
    return None


def extract_upstream_sequence(genome_seqs, contig_name, gene_start, gene_end,
                              strand, upstream_len, head_overlap, logger):
    """
    提取RT基因上游序列（包含头部重叠区）

    坐标换算逻辑：
    - 正链(+): 提取 [Start - upstream, Start + head_overlap]
      序列方向与基因组一致，无需反向互补
    - 负链(-): 提取 [End - head_overlap, End + upstream]
      需要反向互补，使输出为5'->3'编码方向

    返回: (sequence, ext_start, ext_end) 基因组绝对坐标
    """
    contig_record = find_contig(genome_seqs, contig_name)
    if not contig_record:
        logger.debug(f"未找到contig: {contig_name}")
        return None, None, None

    contig_seq = str(contig_record.seq)
    contig_len = len(contig_seq)

    gene_start = max(1, int(gene_start))
    gene_end = min(contig_len, int(gene_end))

    if strand == '+' or strand == '1':
        # 正链：上游在左侧（坐标小的方向）
        # 提取范围: [gene_start - upstream - 1, gene_start + head_overlap - 1]
        ext_start = max(0, gene_start - 1 - upstream_len)
        ext_end = min(contig_len, gene_start - 1 + head_overlap)
        extracted_seq = contig_seq[ext_start:ext_end]
    else:
        # 负链：上游在右侧（坐标大的方向）
        # 提取范围: [gene_end - head_overlap, gene_end + upstream]
        ext_start = max(0, gene_end - head_overlap)
        ext_end = min(contig_len, gene_end + upstream_len)
        extracted_seq = contig_seq[ext_start:ext_end]
        # 负链需要反向互补
        extracted_seq = str(Seq(extracted_seq).reverse_complement())

    return extracted_seq, ext_start, ext_end


def main():
    parser = argparse.ArgumentParser(
        description="Step 2: 上游盲搜序列提取 (Retron ncRNA位于RT基因5'端上游)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python extract_flanking.py -i results.tsv -f /path/to/fasta -o 02_flanking
  python extract_flanking.py -i results.tsv -f /path/to/fasta -o 02_flanking --upstream 800 --head-overlap 100
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Step1输出的TSV文件 (complete_results.tsv)')
    parser.add_argument('-f', '--fasta-dir', required=True,
                        help='基因组fasta文件目录')
    parser.add_argument('-o', '--output', default='02_flanking',
                        help='输出目录 (默认: 02_flanking)')
    parser.add_argument('--upstream', type=int, default=600,
                        help='上游提取长度bp (默认: 600)')
    parser.add_argument('--head-overlap', type=int, default=50,
                        help='RT头部重叠区长度bp (默认: 50)')
    parser.add_argument('--min-identity', type=float, default=30,
                        help='最小identity过滤 (默认: 30)')
    parser.add_argument('--min-coverage', type=float, default=50,
                        help='最小coverage过滤 (默认: 50)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='详细输出')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"输入文件不存在: {input_file}")
        sys.exit(1)

    fasta_dir = Path(args.fasta_dir)
    if not fasta_dir.exists():
        logger.error(f"fasta目录不存在: {fasta_dir}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 2: 上游盲搜序列提取")
    logger.info("=" * 60)
    logger.info(f"提取策略: 上游 {args.upstream}bp + RT头部 {args.head_overlap}bp")
    logger.info("注: Retron ncRNA位于RT基因5'端上游，不提取下游序列")

    # 读取搜索结果
    df = pd.read_csv(input_file, sep='\t')
    logger.info(f"读取到 {len(df)} 条记录")

    if 'identity' in df.columns and 'coverage' in df.columns:
        df = df[(df['identity'] >= args.min_identity) &
                (df['coverage'] >= args.min_coverage)]
        logger.info(f"过滤后: {len(df)} 条记录")

    fasta_index = build_fasta_index(fasta_dir, logger)

    logger.info("\n提取上游序列...")
    records = []
    stats = {'success': 0, 'no_fasta': 0, 'no_contig': 0, 'no_position': 0}
    genome_cache = {}

    for _, row in df.iterrows():
        genome_name = row['genome']

        # 获取位置信息
        if 'genome_nt_start' in df.columns and pd.notna(row.get('genome_nt_start')):
            gene_start = row['genome_nt_start']
            gene_end = row['genome_nt_end']
        elif 'start' in df.columns:
            gene_start = row['start']
            gene_end = row['end']
        else:
            stats['no_position'] += 1
            continue

        # 获取contig名称
        if 'matched_contig' in df.columns and pd.notna(row.get('matched_contig')):
            contig_name = row['matched_contig']
        elif 'contig' in df.columns:
            contig_name = row['contig']
        else:
            stats['no_position'] += 1
            continue

        strand = row.get('gene_strand', row.get('strand', '+'))
        if pd.isna(strand):
            strand = '+'

        fasta_file = find_fasta_file(genome_name, fasta_index)
        if not fasta_file:
            stats['no_fasta'] += 1
            continue

        if str(fasta_file) not in genome_cache:
            genome_cache[str(fasta_file)] = load_genome_sequences(fasta_file)
        genome_seqs = genome_cache[str(fasta_file)]

        seq, ext_start, ext_end = extract_upstream_sequence(
            genome_seqs, contig_name, gene_start, gene_end,
            strand, args.upstream, args.head_overlap, logger
        )

        if seq is None:
            stats['no_contig'] += 1
            continue

        stats['success'] += 1

        # 构建序列ID，包含完整坐标信息用于后续换算
        # 格式: query|genome|protein_id|contig|rt_start-rt_end|strand|ext_start|ext_end
        # protein_id是RT蛋白的唯一标识符（原contig列，如ctg1_5869）
        query_name = row.get('query', 'unknown')
        protein_id = row.get('contig', 'unknown')  # 这是蛋白ID，不是contig
        seq_id = (f"{query_name}|{genome_name}|{protein_id}|{contig_name}|"
                  f"{int(gene_start)}-{int(gene_end)}|{strand}|{ext_start}|{ext_end}")

        record = SeqRecord(
            Seq(seq),
            id=seq_id,
            description=f"upstream_{args.upstream}bp_head_{args.head_overlap}bp"
        )
        records.append(record)

    # 保存结果 - 使用flanking_sequences.fasta保持兼容性
    if records:
        output_file = output_dir / "flanking_sequences.fasta"
        SeqIO.write(records, output_file, "fasta")
        logger.info(f"\n✓ 上游序列: {output_file} ({len(records)} 条)")

    stats_file = output_dir / "extraction_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("上游盲搜序列提取统计\n")
        f.write("=" * 40 + "\n")
        f.write(f"输入记录数: {len(df)}\n")
        f.write(f"成功提取: {stats['success']}\n")
        f.write(f"未找到fasta: {stats['no_fasta']}\n")
        f.write(f"未找到contig: {stats['no_contig']}\n")
        f.write(f"缺少位置信息: {stats['no_position']}\n")
        f.write(f"\n提取参数:\n")
        f.write(f"  上游长度: {args.upstream} bp\n")
        f.write(f"  RT头部重叠: {args.head_overlap} bp\n")
        f.write(f"  总提取长度: {args.upstream + args.head_overlap} bp\n")

    logger.info("\n" + "=" * 60)
    logger.info(f"提取完成! 成功: {stats['success']}/{len(df)}")
    logger.info("=" * 60)

    return 0 if stats['success'] > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
