#!/usr/bin/env python3
"""
Step 3: RT系统发育分析 (Phylogenetic Classification)

基于Mestre et al., 2020的方法，通过系统发育树对RT进行距离标注：
1. 收集Step 2 Motif过滤后的RT序列
2. 添加已知Retron RT参考序列 (Ec67, Ec86, Mx65, Mx162, Stig)
3. MAFFT多序列比对
4. FastTree构建系统发育树
5. 计算到各参考序列的距离，标注置信度

注意：由于没有可靠的Group II/DGR/CRISPR RT参考序列，
本脚本只基于已知Retron RT参考进行距离计算，不进行其他RT类型的分类。

与旧版step03_classify_rt.py的区别：
- 旧版：基于产品描述文本匹配 (不可靠)
- 新版：基于系统发育距离 (科学准确)

输入：Step2的motif_filtered.tsv + FAA文件
输出：
  - retron_candidates.tsv: 所有候选（含距离和置信度标注）
  - rt_classification_summary.tsv: 置信度汇总
  - rt_alignment.fasta: 多序列比对结果
  - rt_tree.nwk: 系统发育树
"""

import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Phylo
import subprocess
import tempfile
import logging
import sys
import re
import os
from collections import defaultdict


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


# 默认参考序列目录
DEFAULT_REFERENCE_DIR = "/home/teng/claude_code/retron/database/rt_reference"


def load_reference_sequences(reference_dir, logger):
    """从目录加载参考RT序列"""
    reference_dir = Path(reference_dir)
    ref_sequences = {}

    if not reference_dir.exists():
        logger.warning(f"参考序列目录不存在: {reference_dir}")
        return ref_sequences

    for fasta_file in reference_dir.glob("*.fasta"):
        try:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                # 使用文件名作为简短ID
                short_name = fasta_file.stem  # e.g., "Ec67_RT_Coli"
                ref_id = f"REF_{short_name}"
                ref_sequences[ref_id] = str(record.seq)
                logger.debug(f"加载参考序列: {ref_id} ({len(record.seq)} aa)")
                break  # 每个文件只取第一条序列
        except Exception as e:
            logger.warning(f"无法解析参考文件 {fasta_file}: {e}")

    logger.info(f"从 {reference_dir} 加载了 {len(ref_sequences)} 条参考序列")
    return ref_sequences


def extract_core_number(name):
    """提取核心数字用于文件匹配"""
    match = re.search(r'(\d{9,}(?:\.\d+)?)', str(name))
    return match.group(1) if match else str(name)


def find_faa_file(genome_name, faa_dir, logger):
    """查找对应的FAA文件"""
    faa_dir = Path(faa_dir)

    # 直接匹配
    for pattern in [f"{genome_name}.faa", f"{genome_name}*.faa"]:
        matches = list(faa_dir.glob(pattern))
        if matches:
            return matches[0]

    # 子目录匹配
    for pattern in [f"**/{genome_name}.faa", f"**/{genome_name}*.faa"]:
        matches = list(faa_dir.rglob(pattern.replace("**/", "")))
        if matches:
            return matches[0]

    # 核心数字匹配
    core = extract_core_number(genome_name)
    for faa_file in faa_dir.rglob("*.faa"):
        if core in str(faa_file):
            return faa_file

    return None


def collect_candidate_sequences(input_file, faa_dir, logger):
    """从Step 2输出收集候选RT序列"""
    df = pd.read_csv(input_file, sep='\t')
    logger.info(f"读取到 {len(df)} 条Motif过滤后的记录")

    sequences = {}
    faa_cache = {}

    for _, row in df.iterrows():
        genome = row.get('genome', '')
        protein_id = row.get('protein_id', row.get('matched_protein', ''))

        if not genome or not protein_id:
            continue

        # 构建唯一ID
        seq_id = f"{genome}|{protein_id}"

        if seq_id in sequences:
            continue

        # 获取FAA文件
        if genome not in faa_cache:
            faa_file = find_faa_file(genome, faa_dir, logger)
            if faa_file:
                faa_cache[genome] = {rec.id: str(rec.seq) for rec in SeqIO.parse(faa_file, 'fasta')}
            else:
                faa_cache[genome] = {}
                logger.debug(f"未找到FAA文件: {genome}")

        # 提取序列
        faa_seqs = faa_cache.get(genome, {})
        seq = None

        # 尝试多种ID格式匹配
        for pid in [protein_id, protein_id.split('|')[0], protein_id.replace('_', ' ')]:
            if pid in faa_seqs:
                seq = faa_seqs[pid]
                break

        # 部分匹配
        if not seq:
            for faa_id, faa_seq in faa_seqs.items():
                if protein_id in faa_id or faa_id in protein_id:
                    seq = faa_seq
                    break

        if seq:
            sequences[seq_id] = {
                'sequence': seq,
                'genome': genome,
                'protein_id': protein_id,
                'row_data': row.to_dict()
            }

    logger.info(f"成功提取 {len(sequences)} 条RT序列")
    return sequences


def add_reference_sequences(sequences, ref_sequences, logger):
    """添加参考RT序列"""
    for ref_id, ref_seq in ref_sequences.items():
        sequences[ref_id] = {
            'sequence': ref_seq,
            'genome': 'REFERENCE',
            'protein_id': ref_id,
            'row_data': {},
            'is_reference': True
        }

    logger.info(f"添加了 {len(ref_sequences)} 条参考序列")
    return sequences


def run_mafft_alignment(sequences, output_dir, logger, threads=8):
    """运行MAFFT多序列比对"""
    input_fasta = output_dir / "rt_sequences_for_alignment.fasta"
    output_fasta = output_dir / "rt_alignment.fasta"

    # 写入序列
    records = []
    for seq_id, seq_data in sequences.items():
        # 清理ID中的特殊字符
        clean_id = seq_id.replace('|', '_').replace(' ', '_')[:80]
        records.append(SeqRecord(
            Seq(seq_data['sequence']),
            id=clean_id,
            description=""
        ))

    SeqIO.write(records, input_fasta, 'fasta')
    logger.info(f"写入 {len(records)} 条序列用于比对")

    # 运行MAFFT
    cmd = [
        'mafft',
        '--auto',
        '--thread', str(threads),
        '--quiet',
        str(input_fasta)
    ]

    logger.info("运行MAFFT比对...")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        with open(output_fasta, 'w') as f:
            f.write(result.stdout)

        logger.info(f"比对完成: {output_fasta}")
        return output_fasta

    except subprocess.CalledProcessError as e:
        logger.error(f"MAFFT运行失败: {e.stderr}")
        return None
    except FileNotFoundError:
        logger.error("MAFFT未安装，请安装: conda install -c bioconda mafft")
        return None


def run_fasttree(alignment_file, output_dir, logger):
    """运行FastTree构建系统发育树"""
    tree_file = output_dir / "rt_tree.nwk"

    cmd = [
        'fasttree',
        '-quiet',
        str(alignment_file)
    ]

    logger.info("运行FastTree构建系统发育树...")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        with open(tree_file, 'w') as f:
            f.write(result.stdout)

        logger.info(f"系统发育树构建完成: {tree_file}")
        return tree_file

    except subprocess.CalledProcessError as e:
        logger.error(f"FastTree运行失败: {e.stderr}")
        return None
    except FileNotFoundError:
        logger.error("FastTree未安装，请安装: conda install -c bioconda fasttree")
        return None


def classify_by_phylogeny(tree_file, sequences, ref_sequences, logger, distance_threshold=0.5):
    """
    基于系统发育树进行RT分类标注

    策略：
    1. 计算每个候选RT到各Retron参考序列的树距离
    2. 找到最近的Retron参考序列
    3. 根据距离远近进行置信度标注

    注意：只使用Retron参考序列进行距离计算，不标注其他RT类型
    """
    try:
        tree = Phylo.read(tree_file, 'newick')
    except Exception as e:
        logger.error(f"无法解析系统发育树: {e}")
        return {}

    # 构建ID映射 (清理后的ID -> 原始ID)
    id_map = {}
    for seq_id in sequences.keys():
        clean_id = seq_id.replace('|', '_').replace(' ', '_')[:80]
        id_map[clean_id] = seq_id

    # 获取所有终端节点
    terminals = {term.name: term for term in tree.get_terminals()}

    # Retron参考序列列表（从动态加载的参考中获取）
    retron_refs = list(ref_sequences.keys())

    all_candidates = {}

    for seq_id, seq_data in sequences.items():
        if seq_data.get('is_reference'):
            continue

        clean_id = seq_id.replace('|', '_').replace(' ', '_')[:80]

        if clean_id not in terminals:
            logger.debug(f"序列不在树中: {seq_id}")
            continue

        query_term = terminals[clean_id]

        # 计算到各Retron参考的距离
        ref_distances = {}
        closest_ref = None
        min_distance = float('inf')

        for ref_id in retron_refs:
            if ref_id in terminals:
                try:
                    dist = tree.distance(query_term, terminals[ref_id])
                    ref_distances[ref_id] = dist
                    if dist < min_distance:
                        min_distance = dist
                        closest_ref = ref_id
                except:
                    pass

        # 基于距离进行置信度标注
        if closest_ref is None:
            # 无法计算距离
            confidence = 'unknown'
            recommended_action = 'review'
        elif min_distance < 0.3:
            # 距离很近，高置信度
            confidence = 'high'
            recommended_action = 'keep'
        elif min_distance < distance_threshold:
            # 距离适中，中等置信度
            confidence = 'medium'
            recommended_action = 'keep'
        elif min_distance < distance_threshold * 1.5:
            # 距离较远但仍在范围内，需要review
            confidence = 'low'
            recommended_action = 'review'
        else:
            # 距离太远，可能不是Retron
            confidence = 'very_low'
            recommended_action = 'review'

        all_candidates[seq_id] = {
            **seq_data,
            'closest_retron_ref': closest_ref if closest_ref else 'N/A',
            'distance_to_closest_ref': min_distance if min_distance != float('inf') else 'N/A',
            'phylo_confidence': confidence,
            'recommended_action': recommended_action,
            'ref_distances': ref_distances,  # 存储所有距离，后续动态输出
        }

    # 统计各置信度数量
    conf_counts = defaultdict(int)
    action_counts = defaultdict(int)
    for data in all_candidates.values():
        conf_counts[data['phylo_confidence']] += 1
        action_counts[data['recommended_action']] += 1

    logger.info(f"系统发育分类完成: 共 {len(all_candidates)} 条")
    logger.info(f"置信度分布:")
    for conf, count in sorted(conf_counts.items()):
        logger.info(f"  {conf}: {count}")
    logger.info(f"建议操作:")
    for action, count in sorted(action_counts.items()):
        logger.info(f"  {action}: {count}")

    return all_candidates


def main():
    parser = argparse.ArgumentParser(
        description="Step 3: RT系统发育分析",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法
  python step03_phylogeny.py -i 02_motif_filtered/motif_filtered.tsv -f /path/to/faa -o 03_phylogeny

  # 指定参考序列目录
  python step03_phylogeny.py -i motif_filtered.tsv -f /path/to/faa -o 03_phylogeny -r /path/to/reference

  # 调整距离阈值
  python step03_phylogeny.py -i motif_filtered.tsv -f /path/to/faa -o 03_phylogeny --distance-threshold 0.6

方法说明:
  基于MAFFT+FastTree构建RT系统发育树，计算候选到已知Retron参考的距离。
  从参考序列目录加载所有.fasta文件作为参考序列。

置信度等级 (phylo_confidence):
  - high: 距离 < 0.3 (高置信度Retron)
  - medium: 0.3 <= 距离 < threshold (中等置信度)
  - low: threshold <= 距离 < threshold*1.5 (低置信度，需review)
  - very_low: 距离 >= threshold*1.5 (可能不是Retron)

建议操作 (recommended_action):
  - keep: 高/中置信度，建议保留
  - review: 低/极低置信度，需要人工检查
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入文件 (Step2的motif_filtered.tsv)')
    parser.add_argument('-f', '--faa-dir', required=True,
                        help='蛋白序列目录 (.faa文件)')
    parser.add_argument('-r', '--reference-dir', default=DEFAULT_REFERENCE_DIR,
                        help=f'参考RT序列目录 (默认: {DEFAULT_REFERENCE_DIR})')
    parser.add_argument('-o', '--output', default='03_phylogeny',
                        help='输出目录 (默认: 03_phylogeny)')
    parser.add_argument('--distance-threshold', type=float, default=0.5,
                        help='系统发育距离阈值，用于判断边界情况 (默认: 0.5)')
    parser.add_argument('--threads', type=int, default=8,
                        help='MAFFT线程数 (默认: 8)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='详细输出')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    # 验证输入
    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"输入文件不存在: {input_file}")
        sys.exit(1)

    faa_dir = Path(args.faa_dir)
    if not faa_dir.exists():
        logger.error(f"FAA目录不存在: {faa_dir}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 3: RT系统发育分析")
    logger.info("=" * 60)
    logger.info(f"输入文件: {input_file}")
    logger.info(f"FAA目录: {faa_dir}")
    logger.info(f"参考序列目录: {args.reference_dir}")
    logger.info(f"输出目录: {output_dir}")

    # 加载参考序列
    logger.info("\n加载参考RT序列...")
    ref_sequences = load_reference_sequences(args.reference_dir, logger)

    if not ref_sequences:
        logger.error("未能加载任何参考序列")
        sys.exit(1)

    # 收集候选序列
    logger.info("\n收集候选RT序列...")
    sequences = collect_candidate_sequences(input_file, faa_dir, logger)

    if not sequences:
        logger.error("未能提取到任何RT序列")
        sys.exit(1)

    # 添加参考序列
    logger.info("\n添加参考RT序列...")
    sequences = add_reference_sequences(sequences, ref_sequences, logger)

    # MAFFT比对
    logger.info("\n运行MAFFT多序列比对...")
    alignment_file = run_mafft_alignment(sequences, output_dir, logger, args.threads)

    if not alignment_file:
        logger.error("MAFFT比对失败")
        sys.exit(1)

    # FastTree建树
    logger.info("\n运行FastTree构建系统发育树...")
    tree_file = run_fasttree(alignment_file, output_dir, logger)

    if not tree_file:
        logger.error("FastTree建树失败")
        sys.exit(1)

    # 系统发育分类
    logger.info("\n基于系统发育树进行分类...")
    all_candidates = classify_by_phylogeny(
        tree_file, sequences, ref_sequences, logger, args.distance_threshold
    )

    # 保存结果
    logger.info("\n保存结果...")

    # 所有候选（含分类标注）
    if all_candidates:
        candidates_data = []
        for seq_id, data in all_candidates.items():
            row = data.get('row_data', {}).copy()
            row['seq_id'] = seq_id
            row['closest_retron_ref'] = data.get('closest_retron_ref', '')
            row['distance_to_closest_ref'] = data.get('distance_to_closest_ref', '')
            row['phylo_confidence'] = data.get('phylo_confidence', '')
            row['recommended_action'] = data.get('recommended_action', '')
            # 动态添加到各参考的距离
            ref_distances = data.get('ref_distances', {})
            for ref_id, dist in ref_distances.items():
                col_name = f"distance_to_{ref_id.replace('REF_', '')}"
                row[col_name] = dist
            candidates_data.append(row)

        candidates_df = pd.DataFrame(candidates_data)
        candidates_file = output_dir / "retron_candidates.tsv"
        candidates_df.to_csv(candidates_file, sep='\t', index=False)
        logger.info(f"  所有候选 (含距离标注): {candidates_file} ({len(candidates_df)} 条)")

        # 置信度汇总表
        summary_data = []
        for conf in candidates_df['phylo_confidence'].unique():
            conf_df = candidates_df[candidates_df['phylo_confidence'] == conf]
            summary_data.append({
                'phylo_confidence': conf,
                'count': len(conf_df),
                'keep': len(conf_df[conf_df['recommended_action'] == 'keep']),
                'review': len(conf_df[conf_df['recommended_action'] == 'review'])
            })

        summary_df = pd.DataFrame(summary_data)
        summary_file = output_dir / "rt_classification_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        logger.info(f"  置信度汇总: {summary_file}")

    # 统计信息
    stats_file = output_dir / "phylogeny_stats.txt"
    distance_threshold = args.distance_threshold  # 获取距离阈值
    with open(stats_file, 'w') as f:
        f.write("RT系统发育分析统计\n")
        f.write("=" * 40 + "\n")
        f.write(f"输入序列数: {len(sequences) - len(ref_sequences)}\n")
        f.write(f"参考序列数: {len(ref_sequences)}\n")
        f.write(f"分类完成数: {len(all_candidates)}\n")
        f.write(f"\n置信度分布:\n")

        if all_candidates:
            conf_counts = defaultdict(int)
            for data in all_candidates.values():
                conf_counts[data.get('phylo_confidence', 'unknown')] += 1
            for conf, count in sorted(conf_counts.items()):
                f.write(f"  {conf}: {count}\n")

            f.write(f"\n建议操作分布:\n")
            action_counts = defaultdict(int)
            for data in all_candidates.values():
                action_counts[data.get('recommended_action', 'unknown')] += 1
            for action, count in sorted(action_counts.items()):
                f.write(f"  {action}: {count}\n")

        f.write(f"\nRetron参考序列 ({len(ref_sequences)}条):\n")
        for ref_id in ref_sequences.keys():
            f.write(f"  - {ref_id}\n")

        f.write(f"\n置信度说明:\n")
        f.write(f"  high: 距离 < 0.3 (高置信度Retron)\n")
        f.write(f"  medium: 0.3 <= 距离 < {distance_threshold} (中等置信度)\n")
        f.write(f"  low: {distance_threshold} <= 距离 < {distance_threshold * 1.5} (低置信度)\n")
        f.write(f"  very_low: 距离 >= {distance_threshold * 1.5} (可能不是Retron)\n")

    # 打印总结
    logger.info("\n" + "=" * 60)
    logger.info("分析完成!")
    logger.info(f"  输入序列: {len(sequences) - len(ref_sequences)}")
    logger.info(f"  分类完成: {len(all_candidates)}")

    if all_candidates:
        conf_counts = defaultdict(int)
        action_counts = defaultdict(int)
        for data in all_candidates.values():
            conf_counts[data.get('phylo_confidence', 'unknown')] += 1
            action_counts[data.get('recommended_action', 'unknown')] += 1

        logger.info(f"\n  置信度分布:")
        for conf, count in sorted(conf_counts.items()):
            logger.info(f"    {conf}: {count}")

        logger.info(f"\n  建议操作:")
        for action, count in sorted(action_counts.items()):
            logger.info(f"    {action}: {count}")

    logger.info(f"\n  输出目录: {output_dir}")
    logger.info("=" * 60)

    return 0 if all_candidates else 1


if __name__ == "__main__":
    sys.exit(main())
