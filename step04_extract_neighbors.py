#!/usr/bin/env python3
"""
Step 6: 全量邻近蛋白提取 (Mestre et al., 2020 方法)

从RT候选上下游 ±30kb 范围内提取所有ORF蛋白序列，
用于后续的MMseqs2聚类和统计关联分析。

与Millman方法的区别：
- Millman: 提取固定数量的邻近基因 (如上下游各5个)
- Mestre: 提取固定距离内的所有ORF (±30kb)，不预设数量限制

输入：Step5筛选结果 + FAA蛋白文件或antiSMASH GBK文件
输出：
  - all_neighbors.faa: 所有邻近蛋白序列 (ID格式: RT_ID|gene_id|position)
  - neighbor_matrix.tsv: RT与邻近基因的对应关系表
  - rt_info.tsv: RT蛋白信息汇总
"""

import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import logging
import sys
import gc
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


def parse_seq_id(seq_id):
    """
    解析序列ID，提取基因组和位置信息

    支持多种格式:
    - 格式1 (Step5): query|genome|protein_id|contig|start-end|strand
    - 格式2 (旧版): query|genome|contig|start-end|strand
    """
    parts = seq_id.split('|')

    # 格式1: 6部分 (Step5 filter_results.tsv)
    if len(parts) >= 6:
        query = parts[0]
        genome = parts[1]
        # parts[2] = protein_id
        contig = parts[3]

        pos_match = re.match(r'(\d+)-(\d+)', parts[4])
        if pos_match:
            start = int(pos_match.group(1))
            end = int(pos_match.group(2))
        else:
            start, end = 0, 0

        strand = parts[5] if len(parts) > 5 else '+'

        return {
            'query': query,
            'genome': genome,
            'contig': contig,
            'start': start,
            'end': end,
            'strand': strand
        }
    # 格式2: 5部分 (旧版)
    elif len(parts) >= 5:
        query = parts[0]
        genome = parts[1]
        contig = parts[2]

        pos_match = re.match(r'(\d+)-(\d+)', parts[3])
        if pos_match:
            start = int(pos_match.group(1))
            end = int(pos_match.group(2))
        else:
            start, end = 0, 0

        strand = parts[4] if len(parts) > 4 else '+'

        return {
            'query': query,
            'genome': genome,
            'contig': contig,
            'start': start,
            'end': end,
            'strand': strand
        }
    return None


def build_antismash_index(antismash_dir, logger):
    """构建antiSMASH目录索引"""
    antismash_dir = Path(antismash_dir)
    index = {}

    for d in antismash_dir.iterdir():
        if d.is_dir():
            index[d.name] = d
            core = extract_core_number(d.name)
            if core not in index:
                index[core] = d

    logger.info(f"索引了 {len(set(index.values()))} 个antiSMASH目录")
    return index


def find_antismash_dir(genome_name, antismash_index):
    """查找对应的antiSMASH目录"""
    if genome_name in antismash_index:
        return antismash_index[genome_name]

    core = extract_core_number(genome_name)
    if core in antismash_index:
        return antismash_index[core]

    return None


# 全局缓存
_gbk_cache = {}


def load_gbk_features(gbk_dir, contig_name, logger):
    """从GBK文件加载CDS特征（带缓存）"""
    global _gbk_cache

    cache_key = (str(gbk_dir), contig_name)
    if cache_key in _gbk_cache:
        return _gbk_cache[cache_key]

    features = []

    gbk_files = [f for f in gbk_dir.glob("*.gbk") if 'region' not in f.name.lower()]

    for gbk_file in gbk_files:
        try:
            for record in SeqIO.parse(gbk_file, 'genbank'):
                record_ids = [record.id, record.name]
                if '.' in record.id:
                    record_ids.append(record.id.split('.')[0])
                if record.id.startswith('NZ_'):
                    record_ids.append(record.id[3:])
                else:
                    record_ids.append('NZ_' + record.id)

                contig_variants = [contig_name]
                if contig_name.startswith('NZ_'):
                    contig_variants.append(contig_name[3:])
                else:
                    contig_variants.append('NZ_' + contig_name)
                if '.' in contig_name:
                    contig_variants.append(contig_name.split('.')[0])

                if not any(cv in record_ids for cv in contig_variants):
                    continue

                for feature in record.features:
                    if feature.type == 'CDS':
                        start = int(feature.location.start) + 1
                        end = int(feature.location.end)
                        strand = '+' if feature.location.strand == 1 else '-'

                        locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                        protein_id = feature.qualifiers.get('protein_id', [''])[0]
                        product = feature.qualifiers.get('product', [''])[0]
                        translation = feature.qualifiers.get('translation', [''])[0]

                        features.append({
                            'locus_tag': locus_tag,
                            'protein_id': protein_id,
                            'product': product,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'translation': translation,
                            'contig': record.id
                        })

        except Exception as e:
            logger.debug(f"解析GBK失败 {gbk_file.name}: {e}")

    features.sort(key=lambda x: x['start'])
    _gbk_cache[cache_key] = features

    return features


def find_neighbors_by_distance(features, target_start, target_end, distance_kb, logger):
    """
    基于距离提取邻近基因 (Mestre方法)

    参数:
        features: 所有CDS特征列表
        target_start, target_end: 目标基因位置
        distance_kb: 上下游提取距离 (kb)

    返回: (upstream_genes, downstream_genes, target_idx)
    """
    if not features:
        return [], [], -1

    distance_bp = distance_kb * 1000
    target_center = (target_start + target_end) / 2

    # 找到目标基因索引
    target_idx = -1
    min_distance = float('inf')

    for i, feat in enumerate(features):
        overlap_start = max(feat['start'], target_start)
        overlap_end = min(feat['end'], target_end)

        if overlap_start <= overlap_end:
            target_idx = i
            break
        else:
            distance = min(abs(feat['start'] - target_end), abs(feat['end'] - target_start))
            if distance < min_distance:
                min_distance = distance
                target_idx = i

    if target_idx < 0:
        return [], [], -1

    # 按距离提取上游基因
    upstream = []
    for i in range(target_idx - 1, -1, -1):
        feat = features[i]
        feat_center = (feat['start'] + feat['end']) / 2
        distance_to_target = target_start - feat['end']

        if distance_to_target > distance_bp:
            break

        upstream.append({
            **feat,
            'distance_to_rt': -distance_to_target,  # 负数表示上游
            'relative_position': 'upstream'
        })

    upstream.reverse()  # 恢复位置顺序

    # 按距离提取下游基因
    downstream = []
    for i in range(target_idx + 1, len(features)):
        feat = features[i]
        distance_to_target = feat['start'] - target_end

        if distance_to_target > distance_bp:
            break

        downstream.append({
            **feat,
            'distance_to_rt': distance_to_target,
            'relative_position': 'downstream'
        })

    return upstream, downstream, target_idx


def create_neighbor_id(rt_id, gene_info):
    """
    创建邻近基因的唯一ID

    格式: RT_ID|locus_tag|position
    这个格式允许后续追踪每个邻近基因属于哪个RT
    """
    locus = gene_info.get('locus_tag', '') or gene_info.get('protein_id', 'unknown')
    pos = gene_info.get('relative_position', 'unknown')
    dist = abs(gene_info.get('distance_to_rt', 0))
    return f"{rt_id}|{locus}|{pos}_{dist}bp"


def main():
    parser = argparse.ArgumentParser(
        description="Step 6: 全量邻近蛋白提取 (Mestre方法)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法
  python step06_extract_neighbors.py -i filter_results.tsv -a /path/to/antismash -o 06_neighbors

  # 自定义提取距离 (默认30kb)
  python step06_extract_neighbors.py -i filter_results.tsv -a /path/to/antismash -o 06_neighbors \\
    --distance 30

  # 从搜索结果直接提取
  python step06_extract_neighbors.py -i complete_results.tsv -a /path/to/antismash -o 06_neighbors \\
    --from-search

方法说明:
  本步骤采用Mestre et al., 2020的方法，提取RT上下游固定距离(±30kb)内的所有ORF，
  而非固定数量的邻近基因。这样可以捕获所有潜在的关联蛋白，用于后续的统计分析。
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入文件 (Step5的filter_results.tsv 或 Step1的complete_results.tsv)')
    parser.add_argument('-a', '--antismash', required=True,
                        help='antiSMASH结果目录')
    parser.add_argument('-o', '--output', default='06_neighbors',
                        help='输出目录 (默认: 06_neighbors)')
    parser.add_argument('--distance', type=int, default=30,
                        help='上下游提取距离，单位kb (默认: 30)')
    parser.add_argument('--from-search', action='store_true',
                        help='输入是Step1的搜索结果（而非Step5的筛选结果）')
    parser.add_argument('--exclude-rt', action='store_true',
                        help='排除RT蛋白本身 (默认包含)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='详细输出')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    # 验证输入
    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"输入文件不存在: {input_file}")
        sys.exit(1)

    antismash_dir = Path(args.antismash)
    if not antismash_dir.exists():
        logger.error(f"antiSMASH目录不存在: {antismash_dir}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 6: 全量邻近蛋白提取 (Mestre方法)")
    logger.info("=" * 60)
    logger.info(f"输入文件: {input_file}")
    logger.info(f"antiSMASH目录: {antismash_dir}")
    logger.info(f"提取距离: ±{args.distance}kb")

    # 构建antiSMASH索引
    logger.info("\n构建antiSMASH目录索引...")
    antismash_index = build_antismash_index(antismash_dir, logger)

    # 读取输入
    logger.info("\n读取输入数据...")
    df = pd.read_csv(input_file, sep='\t')
    logger.info(f"读取到 {len(df)} 条记录")

    # 解析输入格式
    candidates = []

    if args.from_search:
        # Step1 搜索结果格式
        for _, row in df.iterrows():
            candidates.append({
                'query': row.get('query', 'unknown'),
                'genome': row['genome'],
                'contig': row.get('matched_contig', row.get('contig', '')),
                'start': row.get('genome_nt_start', row.get('start', 0)),
                'end': row.get('genome_nt_end', row.get('end', 0)),
                'strand': row.get('gene_strand', row.get('strand', '+'))
            })
    else:
        # Step5 filter_results.tsv 格式 - 优先使用已有的列
        has_direct_columns = all(col in df.columns for col in ['genome', 'contig', 'global_start', 'global_end'])

        if has_direct_columns:
            logger.info("使用独立列: genome, contig, global_start, global_end")
            for _, row in df.iterrows():
                # 从 seq_id 提取 query
                seq_id = row.get('seq_id', '')
                query = seq_id.split('|')[0] if '|' in seq_id else 'unknown'

                candidates.append({
                    'query': query,
                    'genome': row['genome'],
                    'contig': row['contig'],
                    'start': int(row['global_start']),
                    'end': int(row['global_end']),
                    'strand': row.get('strand', '+')
                })
        else:
            # 回退到解析 seq_id
            logger.info("解析 seq_id 列")
            for _, row in df.iterrows():
                seq_id = row['seq_id']
                parsed = parse_seq_id(seq_id)
                if parsed:
                    candidates.append(parsed)

    logger.info(f"解析到 {len(candidates)} 个候选")

    # 按位置去重
    logger.info("\n按位置去重...")
    position_to_candidates = {}
    for cand in candidates:
        pos_key = (cand['genome'], cand['contig'], cand['start'], cand['end'])
        if pos_key not in position_to_candidates:
            position_to_candidates[pos_key] = []
        position_to_candidates[pos_key].append(cand)

    unique_positions = list(position_to_candidates.keys())
    logger.info(f"唯一位置数: {len(unique_positions)} (去重前: {len(candidates)})")

    # 提取邻近蛋白
    logger.info("\n提取邻近蛋白...")

    all_neighbor_proteins = []  # 所有邻近蛋白序列
    neighbor_matrix = []  # RT与邻近基因的对应关系
    rt_info_list = []  # RT信息汇总
    stats = {'success': 0, 'no_dir': 0, 'no_features': 0, 'no_neighbors': 0}

    for idx, pos_key in enumerate(unique_positions):
        genome, contig, target_start, target_end = pos_key
        related_cands = position_to_candidates[pos_key]

        if (idx + 1) % 20 == 0 or idx == 0:
            logger.info(f"处理进度: {idx + 1}/{len(unique_positions)}")

        # 查找antiSMASH目录
        gbk_dir = find_antismash_dir(genome, antismash_index)
        if not gbk_dir:
            stats['no_dir'] += len(related_cands)
            logger.debug(f"未找到antiSMASH目录: {genome}")
            continue

        # 加载GBK特征
        features = load_gbk_features(gbk_dir, contig, logger)
        if not features:
            stats['no_features'] += len(related_cands)
            logger.debug(f"未找到CDS特征: {genome}/{contig}")
            continue

        # 按距离提取邻近基因
        upstream, downstream, target_idx = find_neighbors_by_distance(
            features, target_start, target_end, args.distance, logger
        )

        if target_idx < 0:
            stats['no_neighbors'] += len(related_cands)
            continue

        stats['success'] += len(related_cands)

        # 创建唯一RT ID
        rt_id = f"{genome}_{contig}_{target_start}-{target_end}"
        queries = [c['query'] for c in related_cands]

        # 记录RT信息
        target_gene = features[target_idx] if target_idx < len(features) else None
        rt_info_list.append({
            'rt_id': rt_id,
            'genome': genome,
            'contig': contig,
            'start': target_start,
            'end': target_end,
            'queries': ';'.join(queries),
            'locus_tag': target_gene['locus_tag'] if target_gene else '',
            'product': target_gene['product'] if target_gene else '',
            'upstream_count': len(upstream),
            'downstream_count': len(downstream),
            'total_neighbors': len(upstream) + len(downstream)
        })

        # 收集邻近蛋白
        all_neighbors = upstream + downstream

        # 可选：包含RT本身
        if not args.exclude_rt and target_gene and target_gene.get('translation'):
            target_gene_info = {
                **target_gene,
                'distance_to_rt': 0,
                'relative_position': 'TARGET'
            }
            all_neighbors.append(target_gene_info)

        for gene in all_neighbors:
            if not gene.get('translation'):
                continue

            # 创建唯一的邻近基因ID
            neighbor_id = create_neighbor_id(rt_id, gene)

            # 添加蛋白序列
            all_neighbor_proteins.append(SeqRecord(
                Seq(gene['translation']),
                id=neighbor_id,
                description=f"{gene.get('product', '')} [{gene.get('relative_position', '')} {abs(gene.get('distance_to_rt', 0))}bp]"
            ))

            # 记录对应关系
            neighbor_matrix.append({
                'rt_id': rt_id,
                'neighbor_id': neighbor_id,
                'genome': genome,
                'contig': contig,
                'locus_tag': gene.get('locus_tag', ''),
                'protein_id': gene.get('protein_id', ''),
                'product': gene.get('product', ''),
                'gene_start': gene['start'],
                'gene_end': gene['end'],
                'strand': gene['strand'],
                'distance_to_rt': gene.get('distance_to_rt', 0),
                'relative_position': gene.get('relative_position', '')
            })

        # 定期清理内存
        if (idx + 1) % 50 == 0:
            gc.collect()

    # 保存结果
    logger.info("\n保存结果...")

    # 1. 所有邻近蛋白序列 (用于MMseqs2聚类)
    if all_neighbor_proteins:
        proteins_file = output_dir / "all_neighbors.faa"
        SeqIO.write(all_neighbor_proteins, proteins_file, "fasta")
        logger.info(f"✓ 邻近蛋白序列: {proteins_file} ({len(all_neighbor_proteins)} 条)")

    # 2. 邻近基因对应关系矩阵
    if neighbor_matrix:
        matrix_df = pd.DataFrame(neighbor_matrix)
        matrix_file = output_dir / "neighbor_matrix.tsv"
        matrix_df.to_csv(matrix_file, sep='\t', index=False)
        logger.info(f"✓ 邻近基因矩阵: {matrix_file} ({len(matrix_df)} 条)")

    # 3. RT信息汇总
    if rt_info_list:
        rt_df = pd.DataFrame(rt_info_list)
        rt_file = output_dir / "rt_info.tsv"
        rt_df.to_csv(rt_file, sep='\t', index=False)
        logger.info(f"✓ RT信息汇总: {rt_file} ({len(rt_df)} 个RT)")

    # 4. 统计信息
    stats_file = output_dir / "extraction_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("邻近蛋白提取统计 (Mestre方法)\n")
        f.write("=" * 40 + "\n")
        f.write(f"输入候选数: {len(candidates)}\n")
        f.write(f"唯一位置数: {len(unique_positions)}\n")
        f.write(f"成功提取: {stats['success']}\n")
        f.write(f"未找到antiSMASH目录: {stats['no_dir']}\n")
        f.write(f"未找到CDS特征: {stats['no_features']}\n")
        f.write(f"未找到邻近基因: {stats['no_neighbors']}\n")
        f.write(f"\n输出统计:\n")
        f.write(f"  RT数量: {len(rt_info_list)}\n")
        f.write(f"  邻近蛋白总数: {len(all_neighbor_proteins)}\n")
        f.write(f"  提取距离: ±{args.distance}kb\n")

        if rt_info_list:
            avg_neighbors = sum(r['total_neighbors'] for r in rt_info_list) / len(rt_info_list)
            f.write(f"  平均每个RT的邻近蛋白数: {avg_neighbors:.1f}\n")

    # 打印总结
    logger.info("\n" + "=" * 60)
    logger.info("提取完成!")
    logger.info(f"  唯一RT位置: {len(unique_positions)}")
    logger.info(f"  成功: {stats['success']}/{len(candidates)}")
    logger.info(f"  邻近蛋白总数: {len(all_neighbor_proteins)}")
    if rt_info_list:
        avg_neighbors = sum(r['total_neighbors'] for r in rt_info_list) / len(rt_info_list)
        logger.info(f"  平均每RT邻近蛋白: {avg_neighbors:.1f}")
    logger.info(f"  输出目录: {output_dir}")
    logger.info("=" * 60)

    # 清理缓存
    _gbk_cache.clear()
    gc.collect()

    return 0 if stats['success'] > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
