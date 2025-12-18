#!/usr/bin/env python3
"""
Step 4: 提取邻近基因
从通过RNA筛选的候选中提取RT邻近基因，用于效应蛋白鉴定

输入：Step3筛选结果 + antiSMASH GBK文件
输出：邻近基因蛋白序列fasta + 基因信息表
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

    格式: query|genome|contig|start-end|strand|...
    """
    parts = seq_id.split('|')
    if len(parts) >= 5:
        query = parts[0]
        genome = parts[1]
        contig = parts[2]

        # 解析位置
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


# 全局缓存：{(gbk_dir, contig_name): features}
_gbk_cache = {}

def find_antismash_dir(genome_name, antismash_index):
    """查找对应的antiSMASH目录"""
    if genome_name in antismash_index:
        return antismash_index[genome_name]

    core = extract_core_number(genome_name)
    if core in antismash_index:
        return antismash_index[core]

    return None

def load_gbk_features(gbk_dir, contig_name, logger):
    """从GBK文件加载CDS特征（带缓存）"""
    global _gbk_cache

    # 检查缓存
    cache_key = (str(gbk_dir), contig_name)
    if cache_key in _gbk_cache:
        return _gbk_cache[cache_key]

    features = []

    # 查找GBK文件（排除region文件）
    gbk_files = [f for f in gbk_dir.glob("*.gbk") if 'region' not in f.name.lower()]

    for gbk_file in gbk_files:
        try:
            for record in SeqIO.parse(gbk_file, 'genbank'):
                # 检查contig匹配
                record_ids = [record.id, record.name]
                if '.' in record.id:
                    record_ids.append(record.id.split('.')[0])
                if record.id.startswith('NZ_'):
                    record_ids.append(record.id[3:])
                else:
                    record_ids.append('NZ_' + record.id)

                # 标准化contig_name
                contig_variants = [contig_name]
                if contig_name.startswith('NZ_'):
                    contig_variants.append(contig_name[3:])
                else:
                    contig_variants.append('NZ_' + contig_name)
                if '.' in contig_name:
                    contig_variants.append(contig_name.split('.')[0])

                if not any(cv in record_ids for cv in contig_variants):
                    continue

                # 提取CDS特征
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

    # 按位置排序
    features.sort(key=lambda x: x['start'])

    # 存入缓存
    _gbk_cache[cache_key] = features

    return features

def find_neighbors(features, target_start, target_end, num_neighbors, logger):
    """
    查找目标基因的邻近基因

    返回上游和下游各num_neighbors个基因
    """
    if not features:
        return [], [], -1

    # 找到目标基因的索引
    target_idx = -1
    min_distance = float('inf')

    for i, feat in enumerate(features):
        # 计算与目标位置的重叠或距离
        overlap_start = max(feat['start'], target_start)
        overlap_end = min(feat['end'], target_end)

        if overlap_start <= overlap_end:
            # 有重叠
            target_idx = i
            break
        else:
            # 计算距离
            distance = min(abs(feat['start'] - target_end), abs(feat['end'] - target_start))
            if distance < min_distance:
                min_distance = distance
                target_idx = i

    if target_idx < 0:
        return [], [], -1

    # 提取上游基因
    upstream_start = max(0, target_idx - num_neighbors)
    upstream = features[upstream_start:target_idx]

    # 提取下游基因
    downstream_end = min(len(features), target_idx + num_neighbors + 1)
    downstream = features[target_idx + 1:downstream_end]

    return upstream, downstream, target_idx

def main():
    parser = argparse.ArgumentParser(
        description="Step 4: 提取邻近基因",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 从RNA筛选结果提取邻近基因
  python extract_neighbors.py -i filter_results.tsv -a /path/to/antismash -o 04_neighbors

  # 自定义邻近基因数量
  python extract_neighbors.py -i filter_results.tsv -a /path/to/antismash -o 04_neighbors \\
    --num-neighbors 10

  # 从原始搜索结果提取（跳过RNA筛选）
  python extract_neighbors.py -i complete_results.tsv -a /path/to/antismash -o 04_neighbors \\
    --from-search
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入文件 (Step3的filter_results.tsv 或 Step1的complete_results.tsv)')
    parser.add_argument('-a', '--antismash', required=True,
                        help='antiSMASH结果目录')
    parser.add_argument('-o', '--output', default='04_neighbors',
                        help='输出目录 (默认: 04_neighbors)')
    parser.add_argument('--num-neighbors', type=int, default=5,
                        help='上下游各提取的基因数量 (默认: 5)')
    parser.add_argument('--from-search', action='store_true',
                        help='输入是Step1的搜索结果（而非Step3的筛选结果）')
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

    # 创建输出目录
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 4: 提取邻近基因")
    logger.info("=" * 60)
    logger.info(f"输入文件: {input_file}")
    logger.info(f"antiSMASH目录: {antismash_dir}")
    logger.info(f"邻近基因数量: 上下游各 {args.num_neighbors} 个")

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
        # 从Step1搜索结果读取
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
        # 从Step3筛选结果读取（seq_id格式）
        for _, row in df.iterrows():
            seq_id = row['seq_id']
            parsed = parse_seq_id(seq_id)
            if parsed:
                candidates.append(parsed)

    logger.info(f"解析到 {len(candidates)} 个候选")

    # 去重：按位置分组，避免重复解析
    logger.info("\n按位置去重...")
    position_to_candidates = {}
    for cand in candidates:
        pos_key = (cand['genome'], cand['contig'], cand['start'], cand['end'])
        if pos_key not in position_to_candidates:
            position_to_candidates[pos_key] = []
        position_to_candidates[pos_key].append(cand)

    unique_positions = list(position_to_candidates.keys())
    logger.info(f"唯一位置数: {len(unique_positions)} (去重前: {len(candidates)})")

    # 提取邻近基因
    logger.info("\n提取邻近基因...")

    # 创建clusters目录
    clusters_dir = output_dir / "clusters"
    clusters_dir.mkdir(parents=True, exist_ok=True)

    all_neighbors = []  # 汇总所有邻近基因信息
    cluster_summary = []  # 基因簇汇总
    stats = {'success': 0, 'no_dir': 0, 'no_features': 0, 'no_neighbors': 0}

    for idx, pos_key in enumerate(unique_positions):
        genome, contig, target_start, target_end = pos_key
        related_cands = position_to_candidates[pos_key]

        # 进度报告
        if (idx + 1) % 20 == 0 or idx == 0:
            logger.info(f"处理进度: {idx + 1}/{len(unique_positions)}")

        # 查找antiSMASH目录
        gbk_dir = find_antismash_dir(genome, antismash_index)
        if not gbk_dir:
            stats['no_dir'] += len(related_cands)
            logger.debug(f"未找到antiSMASH目录: {genome}")
            continue

        # 加载GBK特征（使用缓存）
        features = load_gbk_features(gbk_dir, contig, logger)
        if not features:
            stats['no_features'] += len(related_cands)
            logger.debug(f"未找到CDS特征: {genome}/{contig}")
            continue

        # 查找邻近基因
        upstream, downstream, target_idx = find_neighbors(
            features, target_start, target_end, args.num_neighbors, logger
        )

        if target_idx < 0:
            stats['no_neighbors'] += len(related_cands)
            continue

        stats['success'] += len(related_cands)

        # 创建基因簇目录
        cluster_id = f"{genome}_{contig}_{target_start}-{target_end}"
        cluster_dir = clusters_dir / cluster_id
        cluster_dir.mkdir(parents=True, exist_ok=True)

        # 收集该基因簇的邻近基因信息和蛋白序列
        cluster_neighbors = []
        cluster_proteins = []
        target_gene = features[target_idx] if target_idx < len(features) else None

        # 上游基因
        for i, gene in enumerate(upstream):
            position = f"upstream_{len(upstream)-i}"
            cluster_neighbors.append({
                'position': position,
                'distance_to_target': target_start - gene['end'],
                'locus_tag': gene['locus_tag'],
                'protein_id': gene['protein_id'],
                'product': gene['product'],
                'gene_start': gene['start'],
                'gene_end': gene['end'],
                'strand': gene['strand']
            })
            if gene['translation']:
                cluster_proteins.append(SeqRecord(
                    Seq(gene['translation']),
                    id=f"{position}|{gene['locus_tag']}",
                    description=gene['product']
                ))

        # 目标基因 (RT)
        if target_gene:
            cluster_neighbors.append({
                'position': 'TARGET',
                'distance_to_target': 0,
                'locus_tag': target_gene['locus_tag'],
                'protein_id': target_gene['protein_id'],
                'product': target_gene['product'],
                'gene_start': target_gene['start'],
                'gene_end': target_gene['end'],
                'strand': target_gene['strand']
            })
            if target_gene['translation']:
                cluster_proteins.append(SeqRecord(
                    Seq(target_gene['translation']),
                    id=f"TARGET|{target_gene['locus_tag']}",
                    description=target_gene['product']
                ))

        # 下游基因
        for i, gene in enumerate(downstream):
            position = f"downstream_{i+1}"
            cluster_neighbors.append({
                'position': position,
                'distance_to_target': gene['start'] - target_end,
                'locus_tag': gene['locus_tag'],
                'protein_id': gene['protein_id'],
                'product': gene['product'],
                'gene_start': gene['start'],
                'gene_end': gene['end'],
                'strand': gene['strand']
            })
            if gene['translation']:
                cluster_proteins.append(SeqRecord(
                    Seq(gene['translation']),
                    id=f"{position}|{gene['locus_tag']}",
                    description=gene['product']
                ))

        # 保存该基因簇的文件
        # 1. 邻近基因信息
        cluster_info_df = pd.DataFrame(cluster_neighbors)
        cluster_info_df.to_csv(cluster_dir / "neighbors_info.tsv", sep='\t', index=False)

        # 2. 蛋白序列
        if cluster_proteins:
            SeqIO.write(cluster_proteins, cluster_dir / "neighbor_proteins.fasta", "fasta")

        # 3. 元数据
        queries = [c['query'] for c in related_cands]
        with open(cluster_dir / "metadata.txt", 'w') as f:
            f.write(f"Cluster ID: {cluster_id}\n")
            f.write(f"Genome: {genome}\n")
            f.write(f"Contig: {contig}\n")
            f.write(f"Position: {target_start}-{target_end}\n")
            f.write(f"Queries: {', '.join(queries)}\n")
            f.write(f"Upstream genes: {len(upstream)}\n")
            f.write(f"Downstream genes: {len(downstream)}\n")

        # 添加到汇总
        cluster_summary.append({
            'cluster_id': cluster_id,
            'genome': genome,
            'contig': contig,
            'start': target_start,
            'end': target_end,
            'queries': ';'.join(queries),
            'num_queries': len(queries),
            'upstream_count': len(upstream),
            'downstream_count': len(downstream),
            'total_neighbors': len(cluster_neighbors),
            'cluster_dir': str(cluster_dir.relative_to(output_dir))
        })

        # 添加到汇总列表（包含query信息）
        for cand in related_cands:
            for nb in cluster_neighbors:
                record = {
                    'query': cand['query'],
                    'genome': genome,
                    'contig': contig,
                    'cluster_id': cluster_id,
                    **nb
                }
                all_neighbors.append(record)

        # 定期清理内存
        if (idx + 1) % 50 == 0:
            gc.collect()

    # 保存汇总结果
    logger.info("\n保存汇总结果...")

    # 1. 基因簇汇总表
    if cluster_summary:
        summary_df = pd.DataFrame(cluster_summary)
        summary_file = output_dir / "cluster_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        logger.info(f"✓ 基因簇汇总: {summary_file} ({len(summary_df)} 个基因簇)")

    # 2. 所有邻近基因汇总（兼容旧格式）
    if all_neighbors:
        neighbors_df = pd.DataFrame(all_neighbors)
        neighbors_file = output_dir / "neighbors_info.tsv"
        neighbors_df.to_csv(neighbors_file, sep='\t', index=False)
        logger.info(f"✓ 邻近基因汇总: {neighbors_file} ({len(neighbors_df)} 条)")

    # 3. 统计信息
    stats_file = output_dir / "extraction_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("邻近基因提取统计\n")
        f.write("=" * 40 + "\n")
        f.write(f"输入候选数: {len(candidates)}\n")
        f.write(f"唯一位置数: {len(unique_positions)}\n")
        f.write(f"成功提取: {stats['success']}\n")
        f.write(f"未找到antiSMASH目录: {stats['no_dir']}\n")
        f.write(f"未找到CDS特征: {stats['no_features']}\n")
        f.write(f"未找到邻近基因: {stats['no_neighbors']}\n")
        f.write(f"生成基因簇数: {len(cluster_summary)}\n")
        f.write(f"提取的邻近基因总数: {len(all_neighbors)}\n")

    # 打印总结
    logger.info("\n" + "=" * 60)
    logger.info("提取完成!")
    logger.info(f"  唯一位置: {len(unique_positions)}")
    logger.info(f"  生成基因簇: {len(cluster_summary)}")
    logger.info(f"  成功: {stats['success']}/{len(candidates)}")
    logger.info(f"  邻近基因总数: {len(all_neighbors)}")
    logger.info(f"  输出目录: {output_dir}")
    logger.info(f"  基因簇目录: {clusters_dir}")
    logger.info("=" * 60)

    # 清理缓存
    _gbk_cache.clear()
    gc.collect()

    return 0 if stats['success'] > 0 else 1

if __name__ == "__main__":
    sys.exit(main())
