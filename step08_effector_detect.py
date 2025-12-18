#!/usr/bin/env python3
"""
Step 8: 效应蛋白检测

检测Retron候选邻近的效应蛋白特征：
- 2TM (双跨膜) 蛋白
- HNH核酸内切酶
- ATPase
- OLD家族核酸酶
- 核糖基转移酶
- 冷休克蛋白

方法：
1. 关键词匹配（快速，不需要额外工具）
2. InterProScan结构域检测（准确，需要单独环境）
"""

import argparse
import pandas as pd
import subprocess
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


# ============================================================
# Retron相关效应蛋白的InterPro结构域
# ============================================================
EFFECTOR_DOMAINS = {
    # 跨膜蛋白
    'IPR027417': {'type': '2TM_membrane', 'name': 'P-loop NTPase', 'desc': '跨膜相关'},
    'TMhelix': {'type': '2TM_membrane', 'name': 'Transmembrane helix', 'desc': '跨膜螺旋'},

    # 核酸酶
    'IPR003615': {'type': 'HNH_nuclease', 'name': 'HNH nuclease', 'desc': 'HNH核酸内切酶'},
    'IPR002711': {'type': 'HNH_nuclease', 'name': 'HNH endonuclease', 'desc': 'HNH内切酶'},
    'IPR013762': {'type': 'HNH_nuclease', 'name': 'HNH motif', 'desc': 'HNH基序'},
    'IPR004365': {'type': 'OLD_nuclease', 'name': 'OLD family nuclease', 'desc': 'OLD核酸酶'},
    'IPR006309': {'type': 'nuclease', 'name': 'RecB family nuclease', 'desc': 'RecB核酸酶'},

    # ATPase
    'IPR003593': {'type': 'ATPase', 'name': 'AAA+ ATPase', 'desc': 'AAA+ ATPase'},
    'IPR027417': {'type': 'ATPase', 'name': 'P-loop NTPase', 'desc': 'P-loop核苷酸结合'},
    'IPR003439': {'type': 'ATPase', 'name': 'ABC transporter-like', 'desc': 'ABC转运蛋白'},

    # 核糖基转移酶
    'IPR006094': {'type': 'ribosyltransferase', 'name': 'ADP-ribosyltransferase', 'desc': 'ADP核糖基转移酶'},
    'IPR024970': {'type': 'ribosyltransferase', 'name': 'Mono-ADP-ribosyltransferase', 'desc': '单ADP核糖基转移酶'},

    # 冷休克蛋白
    'IPR002059': {'type': 'cold_shock', 'name': 'Cold-shock protein', 'desc': '冷休克蛋白'},
    'IPR011129': {'type': 'cold_shock', 'name': 'CSD domain', 'desc': '冷休克结构域'},

    # 毒素-抗毒素
    'IPR009036': {'type': 'toxin', 'name': 'MazF-like toxin', 'desc': 'MazF毒素'},
    'IPR006442': {'type': 'toxin', 'name': 'RelE/ParE toxin', 'desc': 'RelE/ParE毒素'},
}

# 关键词匹配（备用方案）
EFFECTOR_KEYWORDS = {
    '2TM_membrane': ['transmembrane', 'membrane protein', 'inner membrane', 'integral membrane'],
    'HNH_nuclease': ['HNH', 'endonuclease', 'homing endonuclease', 'DNase'],
    'ATPase': ['ATPase', 'ATP-binding', 'AAA', 'Walker A', 'Walker B'],
    'OLD_nuclease': ['OLD', 'nuclease', 'exonuclease', 'RecB'],
    'ribosyltransferase': ['ribosyltransferase', 'ADP-ribosyl', 'NAD-dependent'],
    'cold_shock': ['cold shock', 'CspA', 'cold-shock'],
    'toxin': ['toxin', 'antitoxin', 'MazF', 'RelE', 'ParE'],
}


def run_interproscan(fasta_file, output_dir, logger, cpu=4):
    """
    运行InterProScan

    返回: 结果TSV文件路径
    """
    output_file = output_dir / "interproscan_results.tsv"

    cmd = [
        'interproscan.sh',
        '-i', str(fasta_file),
        '-f', 'tsv',
        '-o', str(output_file),
        '-cpu', str(cpu),
        '-goterms',
        '-pa',  # Pathway annotations
        '--disable-precalc'  # 不使用预计算结果
    ]

    logger.info(f"运行InterProScan: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200  # 2小时超时
        )

        if result.returncode != 0:
            logger.error(f"InterProScan错误: {result.stderr}")
            return None

        if output_file.exists():
            logger.info(f"InterProScan完成: {output_file}")
            return output_file
        else:
            logger.error("InterProScan未生成输出文件")
            return None

    except subprocess.TimeoutExpired:
        logger.error("InterProScan超时")
        return None
    except FileNotFoundError:
        logger.error("interproscan.sh 未找到，请确保InterProScan已安装并在PATH中")
        return None


def parse_interproscan_results(tsv_file, logger):
    """
    解析InterProScan TSV结果

    TSV格式（11-13列）:
    0: protein_id
    1: md5
    2: length
    3: analysis (Pfam, TIGRFAM, etc.)
    4: signature_accession
    5: signature_description
    6: start
    7: end
    8: score
    9: status
    10: date
    11: interpro_accession (可能为空)
    12: interpro_description (可能为空)
    """
    domain_hits = {}

    try:
        with open(tsv_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue

                protein_id_raw = parts[0]
                # 处理 "upstream_5|ctg1_3054" 格式，提取locus_tag部分
                if '|' in protein_id_raw:
                    protein_id = protein_id_raw.split('|')[-1]
                else:
                    protein_id = protein_id_raw
                analysis = parts[3]
                signature = parts[4]
                sig_desc = parts[5]
                interpro = parts[11] if len(parts) > 11 else ''
                interpro_desc = parts[12] if len(parts) > 12 else ''

                if protein_id not in domain_hits:
                    domain_hits[protein_id] = []

                # 检查是否是效应蛋白相关结构域
                effector_type = None
                effector_info = None

                # 优先检查InterPro ID
                if interpro and interpro in EFFECTOR_DOMAINS:
                    effector_info = EFFECTOR_DOMAINS[interpro]
                    effector_type = effector_info['type']

                # 检查signature
                if not effector_type and signature in EFFECTOR_DOMAINS:
                    effector_info = EFFECTOR_DOMAINS[signature]
                    effector_type = effector_info['type']

                # 检查描述中的关键词
                if not effector_type:
                    desc_lower = (sig_desc + ' ' + interpro_desc).lower()
                    for eff_type, keywords in EFFECTOR_KEYWORDS.items():
                        for kw in keywords:
                            if kw.lower() in desc_lower:
                                effector_type = eff_type
                                effector_info = {'name': signature, 'desc': sig_desc}
                                break
                        if effector_type:
                            break

                if effector_type:
                    domain_hits[protein_id].append({
                        'type': effector_type,
                        'signature': signature,
                        'interpro': interpro,
                        'description': sig_desc or interpro_desc,
                        'analysis': analysis
                    })

        logger.info(f"解析到 {len(domain_hits)} 个蛋白的结构域信息")
        return domain_hits

    except Exception as e:
        logger.error(f"解析InterProScan结果失败: {e}")
        return {}


def classify_by_keywords(product):
    """基于关键词分类（备用方案）"""
    if not product or not isinstance(product, str):
        return []

    product_lower = product.lower()
    matches = []

    for eff_type, keywords in EFFECTOR_KEYWORDS.items():
        for kw in keywords:
            if kw.lower() in product_lower:
                matches.append({
                    'type': eff_type,
                    'method': 'keyword',
                    'match': kw
                })
                break

    return matches


def analyze_effectors(neighbors_df, domain_hits, logger):
    """
    分析邻近基因的效应蛋白特征

    返回: 按候选分组的效应蛋白信息
    """
    effector_info = {}

    for _, row in neighbors_df.iterrows():
        query = row.get('query', 'unknown')
        genome = row['genome']
        position = row.get('position', '')
        product = row.get('product', '')
        locus_tag = row.get('locus_tag', '')
        protein_id = row.get('protein_id', '')

        key = f"{query}|{genome}"

        if key not in effector_info:
            effector_info[key] = {
                'effectors_found': [],
                'effector_types': set(),
                'effector_count': 0,
                'method': 'none'
            }

        # 跳过TARGET基因（RT本身）
        if position == 'TARGET':
            continue

        effector_matches = []

        # 1. 优先使用InterProScan结果
        for pid in [protein_id, locus_tag]:
            if pid and pid in domain_hits:
                for hit in domain_hits[pid]:
                    effector_matches.append({
                        'type': hit['type'],
                        'method': 'interproscan',
                        'signature': hit.get('signature', ''),
                        'description': hit.get('description', '')
                    })
                effector_info[key]['method'] = 'interproscan'

        # 2. 如果没有InterProScan结果，使用关键词匹配
        if not effector_matches:
            keyword_matches = classify_by_keywords(product)
            if keyword_matches:
                effector_matches = keyword_matches
                if effector_info[key]['method'] == 'none':
                    effector_info[key]['method'] = 'keyword'

        # 记录结果
        if effector_matches:
            for match in effector_matches:
                effector_info[key]['effectors_found'].append({
                    'position': position,
                    'locus_tag': locus_tag,
                    'product': product,
                    **match
                })
                effector_info[key]['effector_types'].add(match['type'])

            effector_info[key]['effector_count'] += 1

    return effector_info


def main():
    parser = argparse.ArgumentParser(
        description="Step 8: 效应蛋白检测 (支持InterProScan)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 仅使用关键词匹配（快速）
  python step08_effector_detect.py -i candidates.tsv -n neighbors_info.tsv -o output

  # 使用InterProScan（需要单独环境）
  python step08_effector_detect.py -i candidates.tsv -n neighbors_info.tsv \\
    -p neighbor_proteins.fasta --use-interproscan -o output

  # 使用已有的InterProScan结果
  python step08_effector_detect.py -i candidates.tsv -n neighbors_info.tsv \\
    --interproscan-results interproscan_results.tsv -o output
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='候选结果文件')
    parser.add_argument('-n', '--neighbors', required=True,
                        help='邻近基因信息文件 (neighbors_info.tsv)')
    parser.add_argument('-p', '--proteins',
                        help='邻近基因蛋白序列 (neighbor_proteins.fasta)')
    parser.add_argument('-o', '--output', default='08_effector',
                        help='输出目录')

    # InterProScan选项
    parser.add_argument('--use-interproscan', action='store_true',
                        help='运行InterProScan进行结构域检测')
    parser.add_argument('--interproscan-results',
                        help='使用已有的InterProScan TSV结果文件')
    parser.add_argument('--cpu', type=int, default=4,
                        help='InterProScan使用的CPU数 (默认: 4)')

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    input_file = Path(args.input)
    neighbors_file = Path(args.neighbors)

    if not input_file.exists():
        logger.error(f"输入文件不存在: {input_file}")
        sys.exit(1)
    if not neighbors_file.exists():
        logger.error(f"邻近基因文件不存在: {neighbors_file}")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 8: 效应蛋白检测")
    logger.info("=" * 60)

    # 读取数据
    df = pd.read_csv(input_file, sep='\t')
    neighbors_df = pd.read_csv(neighbors_file, sep='\t')

    logger.info(f"候选数: {len(df)}")
    logger.info(f"邻近基因数: {len(neighbors_df)}")

    # InterProScan分析
    domain_hits = {}

    if args.interproscan_results:
        # 使用已有的结果
        ipr_file = Path(args.interproscan_results)
        if ipr_file.exists():
            logger.info(f"\n使用已有InterProScan结果: {ipr_file}")
            domain_hits = parse_interproscan_results(ipr_file, logger)
        else:
            logger.warning(f"InterProScan结果文件不存在: {ipr_file}")

    elif args.use_interproscan:
        # 运行InterProScan
        if not args.proteins:
            logger.error("使用InterProScan需要提供 -p/--proteins 蛋白序列文件")
            sys.exit(1)

        proteins_file = Path(args.proteins)
        if not proteins_file.exists():
            logger.error(f"蛋白序列文件不存在: {proteins_file}")
            sys.exit(1)

        logger.info(f"\n运行InterProScan...")
        ipr_output = run_interproscan(proteins_file, output_dir, logger, args.cpu)

        if ipr_output:
            domain_hits = parse_interproscan_results(ipr_output, logger)

    if domain_hits:
        logger.info(f"InterProScan检测到 {len(domain_hits)} 个蛋白的结构域")
    else:
        logger.info("\n使用关键词匹配方法")

    # 分析效应蛋白
    logger.info("\n分析效应蛋白...")
    effector_info = analyze_effectors(neighbors_df, domain_hits, logger)

    # 将效应蛋白信息合并到候选
    results = []

    for _, row in df.iterrows():
        query = row.get('query', 'unknown')
        genome = row['genome']
        key = f"{query}|{genome}"

        info = effector_info.get(key, {
            'effector_types': set(),
            'effector_count': 0,
            'effectors_found': [],
            'method': 'none'
        })

        result = row.to_dict()
        result.update({
            'has_effector': info['effector_count'] > 0,
            'effector_count': info['effector_count'],
            'effector_types': ';'.join(info['effector_types']) if info['effector_types'] else '',
            'detection_method': info['method'],
            'effector_details': ';'.join([
                f"{e['position']}:{e['type']}({e.get('method', 'unknown')})"
                for e in info['effectors_found']
            ]) if info['effectors_found'] else ''
        })
        results.append(result)

    result_df = pd.DataFrame(results)

    # 统计
    with_effector = result_df['has_effector'].sum()
    logger.info(f"\n有效应蛋白证据: {with_effector}/{len(result_df)} "
                f"({with_effector/len(result_df)*100:.1f}%)")

    method_counts = result_df['detection_method'].value_counts()
    logger.info(f"检测方法分布: {dict(method_counts)}")

    # 效应蛋白类型统计
    all_types = []
    for types in result_df['effector_types'].dropna():
        if types:
            all_types.extend(types.split(';'))

    if all_types:
        type_counts = pd.Series(all_types).value_counts()
        logger.info("\n效应蛋白类型分布:")
        for eff_type, count in type_counts.items():
            logger.info(f"  {eff_type}: {count}")

    # 保存结果
    output_file = output_dir / "effector_analysis.tsv"
    result_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"\n✓ 结果: {output_file}")

    # 有效应蛋白的候选
    with_effector_df = result_df[result_df['has_effector'] == True]
    if not with_effector_df.empty:
        effector_file = output_dir / "with_effector.tsv"
        with_effector_df.to_csv(effector_file, sep='\t', index=False)
        logger.info(f"✓ 有效应蛋白: {effector_file}")

    # 详细效应蛋白列表
    all_effectors = []
    for key, info in effector_info.items():
        for eff in info['effectors_found']:
            query, genome = key.split('|', 1)
            all_effectors.append({
                'query': query,
                'genome': genome,
                **eff
            })

    if all_effectors:
        effectors_df = pd.DataFrame(all_effectors)
        effectors_file = output_dir / "effector_details.tsv"
        effectors_df.to_csv(effectors_file, sep='\t', index=False)
        logger.info(f"✓ 效应蛋白详情: {effectors_file}")

    logger.info("\n" + "=" * 60)
    logger.info("分析完成!")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
