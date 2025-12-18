#!/usr/bin/env python3
"""
Step 3: 滑动窗口RNA结构筛选 (Sliding Window RNA Filter)

针对黏细菌(High-GC)的Retron挖掘优化版本。
使用滑动窗口扫描 + 严格拓扑筛选，在高GC背景中识别真正的Retron ncRNA。

核心算法：
1. 滑动窗口：对每条序列进行150bp窗口扫描，步长10bp
2. RNAfold：计算每个窗口的二级结构和MFE
3. 拓扑筛选：
   - 颈部闭合：结构两端必须形成配对
   - 分支G：未配对的G必须位于茎结构的气泡中
4. NMS去重：合并重叠窗口，保留MFE最低的

输入：Step2的flanking_sequences.fasta
输出：filter_results.tsv + passed_candidates.fasta
"""

import argparse
import subprocess
import re
import sys
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


def run_rnafold_single(sequence):
    """对单条序列运行RNAfold，返回(structure, mfe)"""
    try:
        result = subprocess.run(
            ['RNAfold', '--noPS'],
            input=f">seq\n{sequence}\n",
            capture_output=True,
            text=True,
            timeout=60
        )
        if result.returncode != 0:
            return None, None

        lines = result.stdout.strip().split('\n')
        if len(lines) >= 3:
            match = re.match(r'^([.()]+)\s*\(([^)]+)\)', lines[2])
            if match:
                return match.group(1), float(match.group(2))
        return None, None
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None, None


def check_neck_closure(structure, min_pairs=3):
    """
    检查结构两端是否形成配对（颈部闭合）

    Retron的msr-msd是独立的结构单元，其根部(a1/a2区域)应该闭合。
    检查前10bp和后10bp的配对数量。
    """
    if len(structure) < 20:
        return False

    left_pairs = structure[:10].count('(')
    right_pairs = structure[-10:].count(')')

    return left_pairs >= min_pairs and right_pairs >= min_pairs


def detect_branching_g_strict(sequence, structure):
    """
    严格版分支G检测：寻找位于茎结构"气泡"中的未配对G

    分支G特征（根据Millman et al., 2020）：
    1. 碱基是G且未配对（结构为'.'）
    2. 位于茎结构的气泡中（前后2bp内有括号）
    3. 高置信度：被茎结构"夹"在中间（左侧有'('且右侧有')'）
    4. 位置限制：不在序列边缘（10%-90%范围内）

    返回: (high_confidence_positions, all_candidate_positions)
    """
    seq_len = len(sequence)
    if seq_len < 20:
        return [], []

    candidates = []
    high_confidence = []

    # 位置限制范围
    min_pos = int(seq_len * 0.1)
    max_pos = int(seq_len * 0.9)

    for i, (base, char) in enumerate(zip(sequence.upper(), structure)):
        if base != 'G' or char != '.':
            continue

        # 位置限制
        if i < min_pos or i > max_pos:
            continue

        # 检查上下文（前后2bp）
        left_context = structure[max(0, i-2):i]
        right_context = structure[i+1:min(seq_len, i+3)]

        left_stem = '(' in left_context
        right_stem = ')' in right_context

        if left_stem or right_stem:
            candidates.append(i)
            # 高置信度：两侧都紧邻茎（气泡结构）
            if left_stem and right_stem:
                high_confidence.append(i)

    return high_confidence, candidates


def find_longest_stem(structure):
    """查找最长的连续茎结构"""
    max_stem = 0
    current = 0
    for char in structure:
        if char == '(':
            current += 1
            max_stem = max(max_stem, current)
        else:
            current = 0
    return max_stem


def sliding_window_scan(sequence, window_size, step_size, min_mfe,
                        min_neck_pairs, require_branch_g, logger):
    """
    滑动窗口扫描序列，返回所有通过筛选的窗口

    返回: [(win_start, win_end, seq, structure, mfe, stem_len, branch_g_count, branch_g_pos), ...]
    """
    seq_len = len(sequence)
    passed_windows = []

    if seq_len < window_size:
        # 序列太短，直接处理整条
        structure, mfe = run_rnafold_single(sequence)
        if structure and mfe and mfe < min_mfe:
            if check_neck_closure(structure, min_neck_pairs):
                high_conf, candidates = detect_branching_g_strict(sequence, structure)
                branch_g_pos = high_conf if high_conf else candidates
                if not require_branch_g or branch_g_pos:
                    stem_len = find_longest_stem(structure)
                    passed_windows.append((0, seq_len, sequence, structure, mfe,
                                          stem_len, len(branch_g_pos), branch_g_pos))
        return passed_windows

    # 滑动窗口扫描
    for start in range(0, seq_len - window_size + 1, step_size):
        end = start + window_size
        window_seq = sequence[start:end]

        structure, mfe = run_rnafold_single(window_seq)
        if not structure or not mfe:
            continue

        # MFE筛选
        if mfe >= min_mfe:
            continue

        # 颈部闭合检测
        if not check_neck_closure(structure, min_neck_pairs):
            continue

        # 分支G检测
        high_conf, candidates = detect_branching_g_strict(window_seq, structure)
        branch_g_pos = high_conf if high_conf else candidates

        if require_branch_g and not branch_g_pos:
            continue

        stem_len = find_longest_stem(structure)
        passed_windows.append((start, end, window_seq, structure, mfe,
                              stem_len, len(branch_g_pos), branch_g_pos))

    return passed_windows


def nms_filter(windows, overlap_threshold=0.5):
    """
    非极大值抑制：合并重叠窗口，保留MFE最低的

    windows: [(start, end, seq, structure, mfe, ...), ...]
    """
    if not windows:
        return []

    # 按MFE排序（从低到高）
    sorted_windows = sorted(windows, key=lambda x: x[4])

    selected = []
    for window in sorted_windows:
        start, end = window[0], window[1]
        window_len = end - start

        # 检查与已选窗口的重叠
        is_overlapping = False
        for sel in selected:
            sel_start, sel_end = sel[0], sel[1]

            # 计算重叠
            overlap_start = max(start, sel_start)
            overlap_end = min(end, sel_end)
            overlap_len = max(0, overlap_end - overlap_start)

            # 重叠比例
            overlap_ratio = overlap_len / min(window_len, sel_end - sel_start)

            if overlap_ratio > overlap_threshold:
                is_overlapping = True
                break

        if not is_overlapping:
            selected.append(window)

    return selected


def calculate_global_coords(win_start, win_end, ext_start, ext_end, strand):
    """
    计算窗口的基因组绝对坐标

    正链: Global = Ext_Start + win_offset
    负链: 序列做了反向互补，需要镜像处理
          Global_Start = Ext_End - win_end
          Global_End = Ext_End - win_start

    返回: (global_start, global_end) 确保 start < end
    """
    if strand == '+' or strand == '1':
        global_start = ext_start + win_start
        global_end = ext_start + win_end
    else:
        # 负链：镜像处理
        global_start = ext_end - win_end
        global_end = ext_end - win_start

    # 确保 start < end
    if global_start > global_end:
        global_start, global_end = global_end, global_start

    return global_start, global_end


def parse_seq_id(seq_id):
    """
    解析序列ID，提取坐标信息
    新格式: query|genome|protein_id|contig|rt_start-rt_end|strand|ext_start|ext_end
    旧格式: query|genome|contig|rt_start-rt_end|strand|ext_start|ext_end
    """
    parts = seq_id.split('|')
    if len(parts) >= 8:
        # 新格式（包含protein_id）
        return {
            'query': parts[0],
            'genome': parts[1],
            'protein_id': parts[2],
            'contig': parts[3],
            'rt_range': parts[4],
            'strand': parts[5],
            'ext_start': int(parts[6]),
            'ext_end': int(parts[7])
        }
    elif len(parts) >= 7:
        # 旧格式（无protein_id）
        return {
            'query': parts[0],
            'genome': parts[1],
            'protein_id': '',
            'contig': parts[2],
            'rt_range': parts[3],
            'strand': parts[4],
            'ext_start': int(parts[5]),
            'ext_end': int(parts[6])
        }
    elif len(parts) >= 5:
        # 更旧的格式
        return {
            'query': parts[0],
            'genome': parts[1],
            'protein_id': '',
            'contig': parts[2],
            'rt_range': parts[3],
            'strand': parts[4] if len(parts) > 4 else '+',
            'ext_start': 0,
            'ext_end': 0
        }
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Step 3: 滑动窗口RNA结构筛选",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python rna_filter.py -i flanking_sequences.fasta -o 03_rna_filtered
  python rna_filter.py -i flanking_sequences.fasta -o 03_rna_filtered --window 150 --step 10 --min-mfe -15
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入fasta文件 (Step2输出)')
    parser.add_argument('-o', '--output', default='03_rna_filtered',
                        help='输出目录 (默认: 03_rna_filtered)')
    parser.add_argument('--window', type=int, default=150,
                        help='滑动窗口大小bp (默认: 150)')
    parser.add_argument('--step', type=int, default=10,
                        help='滑动步长bp (默认: 10)')
    parser.add_argument('--min-mfe', type=float, default=-15.0,
                        help='最小自由能阈值 kcal/mol (默认: -15)')
    parser.add_argument('--min-neck-pairs', type=int, default=3,
                        help='颈部闭合最小配对数 (默认: 3)')
    parser.add_argument('--require-branch-g', action='store_true',
                        help='要求存在分支G')
    # 兼容旧参数
    parser.add_argument('--min-stem', type=int, default=8,
                        help='最小茎长度bp (默认: 8，用于评分)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='详细输出')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"输入文件不存在: {input_file}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 检查RNAfold
    try:
        subprocess.run(['RNAfold', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("RNAfold未安装或不在PATH中")
        logger.error("请安装ViennaRNA: conda install -c bioconda viennarna")
        sys.exit(1)

    logger.info("=" * 60)
    logger.info("Step 3: 滑动窗口RNA结构筛选")
    logger.info("=" * 60)
    logger.info(f"窗口大小: {args.window}bp, 步长: {args.step}bp")
    logger.info(f"MFE阈值: < {args.min_mfe} kcal/mol")
    logger.info(f"颈部闭合: >= {args.min_neck_pairs} 对")
    logger.info(f"要求分支G: {args.require_branch_g}")

    # 读取序列
    sequences = list(SeqIO.parse(input_file, 'fasta'))
    logger.info(f"\n读取到 {len(sequences)} 条序列")

    # 处理每条序列
    all_results = []
    passed_records = []
    stats = {'total': len(sequences), 'passed': 0, 'windows_scanned': 0}

    for i, record in enumerate(sequences):
        if (i + 1) % 10 == 0:
            logger.info(f"处理进度: {i+1}/{len(sequences)}")

        seq_str = str(record.seq)
        parsed = parse_seq_id(record.id)

        if not parsed:
            logger.debug(f"无法解析序列ID: {record.id}")
            continue

        # 滑动窗口扫描
        windows = sliding_window_scan(
            seq_str, args.window, args.step, args.min_mfe,
            args.min_neck_pairs, args.require_branch_g, logger
        )

        stats['windows_scanned'] += max(1, (len(seq_str) - args.window) // args.step + 1)

        if not windows:
            continue

        # NMS去重
        filtered_windows = nms_filter(windows, overlap_threshold=0.5)

        for win in filtered_windows:
            win_start, win_end, win_seq, structure, mfe, stem_len, branch_g_count, branch_g_pos = win

            # 计算基因组绝对坐标
            global_start, global_end = calculate_global_coords(
                win_start, win_end,
                parsed['ext_start'], parsed['ext_end'],
                parsed['strand']
            )

            # 构建兼容后续步骤的seq_id格式
            # 新格式: query|genome|protein_id|contig|global_start-global_end|strand
            protein_id = parsed.get('protein_id', '')
            if protein_id:
                new_seq_id = (f"{parsed['query']}|{parsed['genome']}|{protein_id}|{parsed['contig']}|"
                             f"{global_start}-{global_end}|{parsed['strand']}")
            else:
                # 兼容旧格式
                new_seq_id = (f"{parsed['query']}|{parsed['genome']}|{parsed['contig']}|"
                             f"{global_start}-{global_end}|{parsed['strand']}")

            result = {
                'seq_id': new_seq_id,
                'sequence': win_seq,
                'structure': structure,
                'mfe': mfe,
                'longest_stem': stem_len,
                'branch_g_count': branch_g_count,
                'branch_g_positions': ','.join(map(str, branch_g_pos)) if branch_g_pos else 'none',
                'global_start': global_start,
                'global_end': global_end,
                'genome': parsed['genome'],
                'protein_id': protein_id,
                'contig': parsed['contig'],
                'strand': parsed['strand']
            }
            all_results.append(result)

            # 创建FASTA记录
            passed_record = SeqRecord(
                Seq(win_seq),
                id=new_seq_id,
                description=f"MFE={mfe:.2f} stem={stem_len} branchG={branch_g_count}"
            )
            passed_records.append(passed_record)

        if filtered_windows:
            stats['passed'] += 1

    # 保存结果
    logger.info("\n保存结果...")

    if all_results:
        # 保存TSV（兼容后续步骤）
        results_file = output_dir / "filter_results.tsv"
        with open(results_file, 'w') as f:
            headers = ['seq_id', 'sequence', 'structure', 'mfe', 'longest_stem',
                      'branch_g_count', 'branch_g_positions', 'global_start', 'global_end',
                      'genome', 'protein_id', 'contig', 'strand']
            f.write('\t'.join(headers) + '\n')
            for r in sorted(all_results, key=lambda x: x['mfe']):
                f.write('\t'.join(str(r[h]) for h in headers) + '\n')
        logger.info(f"✓ 筛选结果: {results_file} ({len(all_results)} 条)")

        # 保存FASTA
        passed_fasta = output_dir / "passed_candidates.fasta"
        SeqIO.write(passed_records, passed_fasta, "fasta")
        logger.info(f"✓ 候选序列: {passed_fasta}")

        # 保存结构详情
        struct_file = output_dir / "structures.txt"
        with open(struct_file, 'w') as f:
            for r in sorted(all_results, key=lambda x: x['mfe']):
                f.write(f">{r['seq_id']}\n")
                f.write(f"SEQ: {r['sequence']}\n")
                f.write(f"STR: {r['structure']}\n")
                f.write(f"MFE: {r['mfe']:.2f} kcal/mol\n")
                f.write(f"Stem: {r['longest_stem']} bp\n")
                f.write(f"Branch-G: {r['branch_g_count']} at {r['branch_g_positions']}\n")
                f.write(f"Coords: {r['global_start']}-{r['global_end']} ({r['strand']})\n\n")
        logger.info(f"✓ 结构详情: {struct_file}")

    # 保存统计
    stats_file = output_dir / "filter_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("滑动窗口RNA结构筛选统计\n")
        f.write("=" * 40 + "\n")
        f.write(f"输入序列数: {stats['total']}\n")
        f.write(f"扫描窗口数: {stats['windows_scanned']}\n")
        f.write(f"通过筛选的序列: {stats['passed']}\n")
        f.write(f"通过筛选的窗口: {len(all_results)}\n")
        f.write(f"\n筛选参数:\n")
        f.write(f"  窗口大小: {args.window} bp\n")
        f.write(f"  步长: {args.step} bp\n")
        f.write(f"  MFE阈值: < {args.min_mfe} kcal/mol\n")
        f.write(f"  颈部闭合: >= {args.min_neck_pairs} 对\n")
        f.write(f"  要求分支G: {args.require_branch_g}\n")

    logger.info("\n" + "=" * 60)
    logger.info("筛选完成!")
    logger.info(f"  输入: {stats['total']} 条序列")
    logger.info(f"  通过: {stats['passed']} 条序列, {len(all_results)} 个窗口")
    logger.info(f"  通过率: {stats['passed']/max(1,stats['total'])*100:.1f}%")
    logger.info("=" * 60)

    return 0 if all_results else 1


if __name__ == "__main__":
    sys.exit(main())
