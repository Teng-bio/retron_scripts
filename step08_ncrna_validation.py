#!/usr/bin/env python3
"""
Step 8: ncRNA共变验证 (Covariance Analysis)

基于多序列比对的共变分析验证Retron ncRNA结构：
1. 从候选Retron上游提取核酸序列
2. 按RT系统发育分组（相近的RT应有相似的ncRNA）
3. MAFFT多序列比对
4. RNAalifold共识结构预测
5. R-scape共变分析（统计学验证碱基配对）

与旧版step05_rna_filter.py的区别：
- 旧版：单序列RNAfold（高GC生物产生大量假阳性）
- 新版：多序列共变分析（统计学验证，适用于高GC生物）

输入：Step7候选 + FASTA文件 + Step3系统发育分组
输出：
  - covariance_results.tsv: 共变分析结果
  - group_*/alignment.stk: Stockholm格式比对（含结构）
  - group_*/rscape_results/: R-scape分析结果
  - validated_candidates.tsv: 通过ncRNA验证的候选
"""

import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess
import tempfile
import logging
import sys
import re
import os
import shutil
from collections import defaultdict


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


def find_fasta_file(genome_name, fasta_dir, logger):
    """查找对应的FASTA文件"""
    fasta_dir = Path(fasta_dir)

    # 直接匹配
    for ext in ['.fasta', '.fna', '.fa']:
        for pattern in [f"{genome_name}{ext}", f"{genome_name}*{ext}"]:
            matches = list(fasta_dir.glob(pattern))
            if matches:
                return matches[0]

    # 核心数字匹配
    core = extract_core_number(genome_name)
    for fasta_file in fasta_dir.glob("*"):
        if fasta_file.is_file() and core in str(fasta_file):
            return fasta_file

    return None


def extract_upstream_sequence(genome, contig, start, end, strand, fasta_file,
                              upstream_length=600, logger=None):
    """提取RT上游序列"""
    try:
        # 读取FASTA
        fasta_dict = {}
        for record in SeqIO.parse(fasta_file, 'fasta'):
            fasta_dict[record.id] = record
            # 处理常见ID变体
            if '.' in record.id:
                fasta_dict[record.id.split('.')[0]] = record
            if record.id.startswith('NZ_'):
                fasta_dict[record.id[3:]] = record

        # 查找contig
        contig_record = None
        for cid in [contig, contig.split('.')[0],
                    f"NZ_{contig}" if not contig.startswith('NZ_') else contig[3:]]:
            if cid in fasta_dict:
                contig_record = fasta_dict[cid]
                break

        if not contig_record:
            return None

        contig_seq = str(contig_record.seq)

        # 计算提取坐标
        if strand == '+':
            # 正链：提取start上游
            extract_start = max(0, start - upstream_length - 1)
            extract_end = start - 1
        else:
            # 负链：提取end下游（反向互补后变成上游）
            extract_start = end
            extract_end = min(len(contig_seq), end + upstream_length)

        if extract_start >= extract_end:
            return None

        upstream_seq = contig_seq[extract_start:extract_end]

        # 负链需要反向互补
        if strand == '-':
            upstream_seq = str(Seq(upstream_seq).reverse_complement())

        return upstream_seq

    except Exception as e:
        if logger:
            logger.debug(f"提取序列失败: {e}")
        return None


def group_candidates_by_phylogeny(candidates_df, phylo_file, logger):
    """
    基于系统发育分组候选

    如果没有系统发育信息，则基于query类型分组
    """
    groups = defaultdict(list)

    if phylo_file and Path(phylo_file).exists():
        # 读取系统发育分类结果
        phylo_df = pd.read_csv(phylo_file, sep='\t')

        # 基于closest_ref分组
        for _, row in candidates_df.iterrows():
            seq_id = row.get('seq_id', f"{row.get('genome', '')}|{row.get('protein_id', '')}")

            # 查找对应的系统发育分类
            closest_ref = 'unknown'
            for _, pr in phylo_df.iterrows():
                if pr.get('seq_id', '') == seq_id or \
                   (pr.get('genome', '') == row.get('genome', '') and
                    pr.get('protein_id', '') == row.get('protein_id', '')):
                    closest_ref = str(pr.get('closest_ref', 'unknown'))
                    break

            groups[closest_ref].append(row.to_dict())
    else:
        # 没有系统发育文件，基于query分组
        for _, row in candidates_df.iterrows():
            query = row.get('query', 'unknown')
            groups[query].append(row.to_dict())

    logger.info(f"候选分为 {len(groups)} 组")
    for group_name, members in groups.items():
        logger.debug(f"  {group_name}: {len(members)} 条")

    return dict(groups)


def run_mafft_alignment(sequences, output_file, logger, threads=4):
    """运行MAFFT多序列比对"""
    if len(sequences) < 2:
        logger.debug("序列数少于2，跳过比对")
        return None

    # 写入临时文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        for seq_id, seq in sequences.items():
            tmp.write(f">{seq_id}\n{seq}\n")
        tmp_file = tmp.name

    try:
        cmd = [
            'mafft',
            '--auto',
            '--thread', str(threads),
            '--quiet',
            tmp_file
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        with open(output_file, 'w') as f:
            f.write(result.stdout)

        return output_file

    except Exception as e:
        logger.debug(f"MAFFT失败: {e}")
        return None
    finally:
        os.unlink(tmp_file)


def run_rnaalifold(alignment_file, output_dir, logger):
    """运行RNAalifold获取共识结构"""
    try:
        # RNAalifold需要Stockholm格式，先转换
        stk_file = output_dir / "alignment.stk"

        # 读取FASTA格式比对，转换为Stockholm
        records = list(SeqIO.parse(alignment_file, 'fasta'))
        if len(records) < 2:
            return None

        # 写入Stockholm格式
        with open(stk_file, 'w') as f:
            f.write("# STOCKHOLM 1.0\n\n")
            for rec in records:
                # Stockholm格式: name sequence
                f.write(f"{rec.id[:30]:<30} {str(rec.seq)}\n")
            f.write("//\n")

        # 运行RNAalifold
        cmd = [
            'RNAalifold',
            '--noPS',
            '-q',
            str(stk_file)
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(output_dir)
        )

        # 解析输出
        lines = result.stdout.strip().split('\n')
        if len(lines) >= 2:
            consensus_seq = lines[0].strip()
            structure_line = lines[1].strip()

            # 解析结构和能量
            # 格式: "(((...))) (-10.50)"
            parts = structure_line.rsplit(' ', 1)
            if len(parts) == 2:
                structure = parts[0]
                mfe_match = re.search(r'\(([-\d.]+)\)', parts[1])
                mfe = float(mfe_match.group(1)) if mfe_match else 0.0

                return {
                    'consensus_seq': consensus_seq,
                    'structure': structure,
                    'mfe': mfe,
                    'stk_file': str(stk_file)
                }

        return None

    except Exception as e:
        logger.debug(f"RNAalifold失败: {e}")
        return None


def run_rscape(stk_file, output_dir, logger):
    """
    运行R-scape进行共变分析

    R-scape检测RNA多序列比对中的统计显著碱基配对
    """
    try:
        cmd = [
            'R-scape',
            '--outdir', str(output_dir),
            '-s',  # 输出统计信息
            str(stk_file)
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(output_dir)
        )

        # 解析R-scape输出
        rscape_results = {
            'significant_pairs': 0,
            'e_value': None,
            'covariation_score': 0,
            'output_files': []
        }

        # 查找输出文件
        for f in output_dir.glob("*.cov"):
            rscape_results['output_files'].append(str(f))

            # 解析.cov文件统计显著配对
            with open(f) as cov_file:
                for line in cov_file:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 6:
                        # 格式: i j score E-value ...
                        try:
                            e_val = float(parts[3])
                            if e_val < 0.05:  # 显著性阈值
                                rscape_results['significant_pairs'] += 1
                        except:
                            pass

        # 从stdout解析额外信息
        for line in result.stdout.split('\n'):
            if 'E-value' in line:
                match = re.search(r'E-value:\s*([\d.e-]+)', line)
                if match:
                    rscape_results['e_value'] = float(match.group(1))

        return rscape_results

    except FileNotFoundError:
        logger.warning("R-scape未安装，跳过共变分析")
        logger.warning("安装: conda install -c bioconda r-scape")
        return None
    except Exception as e:
        logger.debug(f"R-scape失败: {e}")
        return None


def analyze_group(group_name, members, fasta_dir, output_dir, logger, upstream_length=600):
    """分析单个分组的ncRNA"""
    group_dir = output_dir / f"group_{group_name.replace('/', '_').replace(' ', '_')[:30]}"
    group_dir.mkdir(parents=True, exist_ok=True)

    # 提取上游序列
    sequences = {}
    fasta_cache = {}

    for member in members:
        genome = member.get('genome', '')
        contig = member.get('contig', member.get('matched_contig', ''))
        start = member.get('global_start', member.get('start', 0))
        end = member.get('global_end', member.get('end', 0))
        strand = member.get('strand', '+')

        if not genome or not contig:
            continue

        # 获取FASTA文件
        if genome not in fasta_cache:
            fasta_file = find_fasta_file(genome, fasta_dir, logger)
            fasta_cache[genome] = fasta_file

        fasta_file = fasta_cache.get(genome)
        if not fasta_file:
            continue

        # 提取序列
        seq = extract_upstream_sequence(
            genome, contig, int(start), int(end), strand,
            fasta_file, upstream_length, logger
        )

        if seq and len(seq) >= 100:  # 最小长度要求
            seq_id = f"{genome}_{contig}_{start}-{end}"
            sequences[seq_id] = seq

    if len(sequences) < 3:
        logger.debug(f"组 {group_name}: 序列数不足 ({len(sequences)})")
        return {
            'group': group_name,
            'n_sequences': len(sequences),
            'status': 'insufficient_sequences'
        }

    # MAFFT比对
    alignment_file = group_dir / "alignment.fasta"
    aligned = run_mafft_alignment(sequences, alignment_file, logger)

    if not aligned:
        return {
            'group': group_name,
            'n_sequences': len(sequences),
            'status': 'alignment_failed'
        }

    # RNAalifold共识结构
    alifold_result = run_rnaalifold(alignment_file, group_dir, logger)

    if not alifold_result:
        return {
            'group': group_name,
            'n_sequences': len(sequences),
            'status': 'alifold_failed'
        }

    # R-scape共变分析
    rscape_result = run_rscape(alifold_result['stk_file'], group_dir, logger)

    result = {
        'group': group_name,
        'n_sequences': len(sequences),
        'status': 'success',
        'consensus_structure': alifold_result['structure'],
        'mfe': alifold_result['mfe'],
        'alignment_file': str(alignment_file),
        'members': [m.get('seq_id', f"{m.get('genome')}|{m.get('protein_id')}") for m in members]
    }

    if rscape_result:
        result['significant_pairs'] = rscape_result['significant_pairs']
        result['covariation_validated'] = rscape_result['significant_pairs'] >= 2

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Step 8: ncRNA共变验证",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法
  python step08_ncrna_validation.py -i 07_report/final_candidates.tsv -f /path/to/fasta -o 08_ncrna

  # 指定系统发育分组
  python step08_ncrna_validation.py -i final_candidates.tsv -f /path/to/fasta -o 08_ncrna \\
    --phylogeny 03_phylogeny/retron_candidates.tsv

方法说明:
  基于多序列比对的共变分析，使用RNAalifold和R-scape验证ncRNA结构。
  这种方法对高GC生物更可靠，因为它依赖统计共变而非单序列折叠。
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入文件 (Step7的final_candidates.tsv)')
    parser.add_argument('-f', '--fasta-dir', required=True,
                        help='基因组核酸序列目录 (.fasta文件)')
    parser.add_argument('-o', '--output', default='08_ncrna',
                        help='输出目录 (默认: 08_ncrna)')
    parser.add_argument('--phylogeny',
                        help='系统发育分类结果 (Step3的retron_candidates.tsv)')
    parser.add_argument('--upstream', type=int, default=600,
                        help='上游提取长度 (默认: 600bp)')
    parser.add_argument('--min-group-size', type=int, default=3,
                        help='最小分组大小 (默认: 3)')
    parser.add_argument('--threads', type=int, default=4,
                        help='MAFFT线程数 (默认: 4)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='详细输出')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    # 验证输入
    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"输入文件不存在: {input_file}")
        sys.exit(1)

    fasta_dir = Path(args.fasta_dir)
    if not fasta_dir.exists():
        logger.error(f"FASTA目录不存在: {fasta_dir}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 8: ncRNA共变验证")
    logger.info("=" * 60)
    logger.info(f"输入文件: {input_file}")
    logger.info(f"FASTA目录: {fasta_dir}")
    logger.info(f"输出目录: {output_dir}")

    # 检查工具
    tools_available = {
        'mafft': shutil.which('mafft') is not None,
        'RNAalifold': shutil.which('RNAalifold') is not None,
        'R-scape': shutil.which('R-scape') is not None
    }

    for tool, available in tools_available.items():
        if available:
            logger.info(f"  {tool}: 已安装")
        else:
            logger.warning(f"  {tool}: 未安装")

    if not tools_available['mafft']:
        logger.error("MAFFT是必需的，请安装: conda install -c bioconda mafft")
        sys.exit(1)

    # 读取候选
    df = pd.read_csv(input_file, sep='\t')
    logger.info(f"\n读取到 {len(df)} 条候选")

    # 分组
    logger.info("\n按系统发育分组...")
    groups = group_candidates_by_phylogeny(df, args.phylogeny, logger)

    # 分析每个组
    logger.info("\n分析各组ncRNA...")
    all_results = []

    for group_name, members in groups.items():
        if len(members) < args.min_group_size:
            logger.debug(f"跳过组 {group_name}: 成员数 {len(members)} < {args.min_group_size}")
            continue

        logger.info(f"\n处理组: {group_name} ({len(members)} 条)")
        result = analyze_group(
            group_name, members, fasta_dir, output_dir,
            logger, args.upstream
        )
        all_results.append(result)

        if result['status'] == 'success':
            sig_pairs = result.get('significant_pairs', 0)
            logger.info(f"  共识结构MFE: {result.get('mfe', 'N/A')}")
            logger.info(f"  显著共变配对: {sig_pairs}")

    # 汇总结果
    logger.info("\n保存结果...")

    if all_results:
        # 共变分析结果
        results_df = pd.DataFrame(all_results)
        results_file = output_dir / "covariance_results.tsv"
        results_df.to_csv(results_file, sep='\t', index=False)
        logger.info(f"  共变分析结果: {results_file}")

        # 通过验证的候选
        validated_groups = [r for r in all_results
                           if r.get('covariation_validated', False)]

        if validated_groups:
            # 提取通过验证的成员
            validated_ids = set()
            for vg in validated_groups:
                validated_ids.update(vg.get('members', []))

            # 筛选原始数据
            validated_candidates = []
            for _, row in df.iterrows():
                seq_id = row.get('seq_id', f"{row.get('genome')}|{row.get('protein_id')}")
                if seq_id in validated_ids:
                    row_dict = row.to_dict()
                    row_dict['ncrna_validated'] = True
                    validated_candidates.append(row_dict)

            if validated_candidates:
                validated_df = pd.DataFrame(validated_candidates)
                validated_file = output_dir / "validated_candidates.tsv"
                validated_df.to_csv(validated_file, sep='\t', index=False)
                logger.info(f"  验证通过的候选: {validated_file} ({len(validated_df)} 条)")

    # 统计
    stats_file = output_dir / "validation_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("ncRNA共变验证统计\n")
        f.write("=" * 40 + "\n")
        f.write(f"输入候选数: {len(df)}\n")
        f.write(f"分组数: {len(groups)}\n")
        f.write(f"分析的组数: {len(all_results)}\n")

        success_count = sum(1 for r in all_results if r['status'] == 'success')
        validated_count = sum(1 for r in all_results if r.get('covariation_validated', False))

        f.write(f"成功分析: {success_count}\n")
        f.write(f"通过共变验证: {validated_count}\n")

        if all_results:
            f.write(f"\n各组详情:\n")
            for r in all_results:
                f.write(f"  {r['group']}: {r['status']}")
                if r['status'] == 'success':
                    f.write(f", MFE={r.get('mfe', 'N/A')}, sig_pairs={r.get('significant_pairs', 0)}")
                f.write("\n")

    # 打印总结
    logger.info("\n" + "=" * 60)
    logger.info("验证完成!")
    logger.info(f"  输入候选: {len(df)}")
    logger.info(f"  分析的组: {len(all_results)}")

    success_count = sum(1 for r in all_results if r['status'] == 'success')
    validated_count = sum(1 for r in all_results if r.get('covariation_validated', False))

    logger.info(f"  成功分析: {success_count}")
    logger.info(f"  通过共变验证: {validated_count}")
    logger.info(f"  输出目录: {output_dir}")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
