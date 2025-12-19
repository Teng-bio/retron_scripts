#!/usr/bin/env python3
"""
Step 1 (HMMER版本): 使用自建Retron RT HMM模型搜索候选序列

与step01_search.py类似的流式处理架构，但使用hmmsearch代替DIAMOND/BLAST
优势：基于Retron RT特异的HMM模型，能准确识别RT核心区域边界(env_from/env_to)

输入：
  - HMM模型文件 (Retron_RT.hmm)
  - FAA蛋白序列目录

输出：
  - complete_results.tsv (与step01格式兼容)
  - 包含 env_from/env_to 列用于后续核心区域提取
"""

import os
import sys
import subprocess
import logging
import shutil
from pathlib import Path
import time
from Bio import SeqIO
import pandas as pd
import argparse
import gc
import psutil
import re
from collections import defaultdict
import tempfile


def get_memory_usage():
    """获取当前内存使用情况（MB）"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024


def setup_logging(output_dir, verbose=False):
    """设置日志系统"""
    log_file = output_dir / "logs" / "hmmsearch_pipeline.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)

    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def extract_contig_from_protein_id(protein_id):
    """从protein_id中提取contig名称"""
    protein_id = protein_id.split()[0]

    if '|' in protein_id:
        parts = protein_id.rsplit('_', 1)
        if len(parts) > 1 and parts[-1].isdigit():
            return parts[0]

    if '_' in protein_id:
        parts = protein_id.rsplit('_', 1)
        if parts[-1].isdigit():
            return parts[0]

    match = re.match(r'([A-Z_]+\d+\.\d+)', protein_id)
    if match:
        return match.group(1)

    return protein_id


def parse_sequence_location(description):
    """从FASTA描述中解析位置信息"""
    patterns = [
        r'#\s*(\d+)\s*#\s*(\d+)\s*#',
        r'Location[:\s]+([\d,]+)\s*-\s*([\d,]+)',
        r'location[=:]\s*(\d+)\.\.(\d+)',
    ]

    for pattern in patterns:
        match = re.search(pattern, description, re.IGNORECASE)
        if match:
            start = int(match.group(1).replace(',', ''))
            end = int(match.group(2).replace(',', ''))
            return start, end

    return None, None


def parse_strand_from_description(description):
    """从描述中解析链方向"""
    match = re.search(r'#\s*\d+\s*#\s*\d+\s*#\s*([-+]?\d+|[-+])', description)
    if match:
        strand_str = match.group(1)
        if strand_str in ['1', '+']:
            return '+'
        elif strand_str in ['-1', '-']:
            return '-'
    return None


def run_hmmsearch(hmm_file, fasta_file, output_dir, threads=4, evalue=1e-5):
    """运行hmmsearch并解析结果"""
    domtblout = output_dir / "hmmsearch_domtblout.txt"
    tblout = output_dir / "hmmsearch_tblout.txt"
    full_output = output_dir / "hmmsearch.out"

    cmd = [
        'hmmsearch',
        '--domtblout', str(domtblout),
        '--tblout', str(tblout),
        '-o', str(full_output),
        '--cpu', str(threads),
        '-E', str(evalue),
        str(hmm_file),
        str(fasta_file)
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=3600)
        return domtblout
    except subprocess.CalledProcessError as e:
        logging.error(f"hmmsearch失败: {e.stderr.decode()}")
        return None
    except FileNotFoundError:
        logging.error("hmmsearch未安装，请安装: conda install -c bioconda hmmer")
        return None


def parse_domtblout(domtblout_file):
    """解析hmmsearch domtblout输出"""
    hits = []

    with open(domtblout_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 23:
                continue

            try:
                hit = {
                    'target_name': parts[0],
                    'target_accession': parts[1],
                    'target_len': int(parts[2]),
                    'query_name': parts[3],
                    'query_accession': parts[4],
                    'query_len': int(parts[5]),
                    'full_evalue': float(parts[6]),
                    'full_score': float(parts[7]),
                    'full_bias': float(parts[8]),
                    'domain_num': int(parts[9]),
                    'domain_total': int(parts[10]),
                    'domain_cevalue': float(parts[11]),
                    'domain_ievalue': float(parts[12]),
                    'domain_score': float(parts[13]),
                    'domain_bias': float(parts[14]),
                    'hmm_from': int(parts[15]),
                    'hmm_to': int(parts[16]),
                    'ali_from': int(parts[17]),
                    'ali_to': int(parts[18]),
                    'env_from': int(parts[19]),
                    'env_to': int(parts[20]),
                    'accuracy': float(parts[21]),
                    'description': ' '.join(parts[22:]) if len(parts) > 22 else ''
                }
                hits.append(hit)
            except (ValueError, IndexError) as e:
                continue

    return hits


class HMMSearchPipeline:
    """HMMER搜索流程"""

    def __init__(self, hmm_file, database_dir, output_dir="hmm_search_results",
                 batch_size=100, threads=4, evalue=1e-5, min_score=25.0):

        self.hmm_file = Path(hmm_file)
        self.database_dir = Path(database_dir)
        self.output_dir = Path(output_dir)
        self.batch_size = batch_size
        self.threads = threads
        self.evalue = evalue
        self.min_score = min_score

        self.temp_dir = Path(tempfile.mkdtemp(prefix="hmmsearch_"))

        self._create_output_structure()
        self.logger = setup_logging(self.output_dir)
        self._validate_inputs()

        self.database_info = self._scan_database_files()
        self.batch_plan = self._create_batch_plan()

    def _create_output_structure(self):
        """创建输出目录结构"""
        directories = ["reports", "logs", "batches", "hmmsearch"]
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for dir_name in directories:
            (self.output_dir / dir_name).mkdir(parents=True, exist_ok=True)

    def _validate_inputs(self):
        """验证输入文件"""
        if not self.hmm_file.exists():
            self.logger.error(f"HMM模型文件不存在: {self.hmm_file}")
            sys.exit(1)

        if not self.database_dir.exists():
            self.logger.error(f"数据库目录不存在: {self.database_dir}")
            sys.exit(1)

        if not shutil.which('hmmsearch'):
            self.logger.error("hmmsearch未安装，请安装: conda install -c bioconda hmmer")
            sys.exit(1)

        self.logger.info(f"HMM模型: {self.hmm_file}")
        self.logger.info(f"数据库目录: {self.database_dir}")

    def _scan_database_files(self):
        """扫描数据库文件"""
        self.logger.info("扫描数据库文件...")

        fasta_extensions = ['*.faa', '*.fasta', '*.fa', '*.pep']
        db_files = []

        for ext in fasta_extensions:
            db_files.extend(self.database_dir.rglob(ext))

        if not db_files:
            self.logger.error(f"未找到序列文件: {self.database_dir}")
            sys.exit(1)

        database_info = {}
        for db_file in db_files:
            genome_name = db_file.stem
            database_info[genome_name] = {'file': db_file}

        self.logger.info(f"找到 {len(database_info)} 个序列文件")
        return database_info

    def _create_batch_plan(self):
        """创建批次计划"""
        db_names = list(self.database_info.keys())

        batches = []
        for i in range(0, len(db_names), self.batch_size):
            batch_items = db_names[i:i + self.batch_size]
            batches.append({
                'batch_id': i // self.batch_size + 1,
                'items': batch_items,
                'size': len(batch_items)
            })

        self.logger.info(f"创建 {len(batches)} 个批次")
        return batches

    def _create_batch_fasta(self, batch_items, batch_id):
        """创建批次FASTA文件"""
        batch_fasta = self.temp_dir / f"batch_{batch_id}.faa"
        genome_mapping = {}
        positions = {}

        with open(batch_fasta, 'w') as out_f:
            for genome_name in batch_items:
                info = self.database_info[genome_name]
                fasta_file = info['file']

                for record in SeqIO.parse(fasta_file, 'fasta'):
                    original_id = record.id
                    new_id = f"{genome_name}___{original_id}"

                    # 提取位置信息
                    nt_start, nt_end = parse_sequence_location(record.description)
                    strand = parse_strand_from_description(record.description)
                    contig = extract_contig_from_protein_id(original_id)

                    positions[new_id] = {
                        'nt_start': nt_start,
                        'nt_end': nt_end,
                        'strand': strand,
                        'contig': contig,
                        'aa_length': len(record.seq),
                        'description': record.description
                    }

                    genome_mapping[new_id] = {
                        'genome': genome_name,
                        'sequence_id': original_id
                    }

                    out_f.write(f">{new_id}\n{record.seq}\n")

        return batch_fasta, genome_mapping, positions

    def _process_batch(self, batch):
        """处理单个批次"""
        batch_id = batch['batch_id']
        batch_items = batch['items']

        self.logger.info(f"\n{'='*20} 批次 {batch_id}/{len(self.batch_plan)} {'='*20}")
        self.logger.info(f"包含 {len(batch_items)} 个基因组")

        # 创建批次FASTA
        batch_fasta, genome_mapping, positions = self._create_batch_fasta(batch_items, batch_id)

        # 运行hmmsearch
        self.logger.info("运行 hmmsearch...")
        batch_output = self.output_dir / "hmmsearch" / f"batch_{batch_id}"
        batch_output.mkdir(parents=True, exist_ok=True)

        domtblout = run_hmmsearch(
            self.hmm_file, batch_fasta, batch_output,
            self.threads, self.evalue
        )

        if not domtblout or not domtblout.exists():
            self.logger.warning(f"批次 {batch_id} hmmsearch失败")
            return None

        # 解析结果
        hits = parse_domtblout(domtblout)
        self.logger.info(f"找到 {len(hits)} 个domain匹配")

        if not hits:
            return None

        # 转换为DataFrame
        results = []
        for hit in hits:
            target = hit['target_name']

            if target not in genome_mapping:
                continue

            info = genome_mapping[target]
            pos_info = positions.get(target, {})

            # 过滤低分匹配
            if hit['domain_score'] < self.min_score:
                continue

            result = {
                'query': 'Retron_RT_HMM',
                'genome': info['genome'],
                'protein_id': info['sequence_id'],
                'contig': pos_info.get('contig', ''),

                # HMM匹配信息
                'hmm_score': hit['domain_score'],
                'hmm_evalue': hit['domain_ievalue'],
                'hmm_bias': hit['domain_bias'],
                'accuracy': hit['accuracy'],

                # 核心区域边界 (关键!)
                'env_from': hit['env_from'],
                'env_to': hit['env_to'],
                'env_length': hit['env_to'] - hit['env_from'] + 1,

                # 比对区域
                'ali_from': hit['ali_from'],
                'ali_to': hit['ali_to'],

                # HMM覆盖
                'hmm_from': hit['hmm_from'],
                'hmm_to': hit['hmm_to'],
                'hmm_coverage': (hit['hmm_to'] - hit['hmm_from'] + 1) / hit['query_len'] * 100,

                # 序列信息
                'target_len': hit['target_len'],
                'domain_num': hit['domain_num'],
                'domain_total': hit['domain_total'],

                # 位置信息
                'nt_start': pos_info.get('nt_start'),
                'nt_end': pos_info.get('nt_end'),
                'strand': pos_info.get('strand'),

                'method': 'HMMSEARCH'
            }
            results.append(result)

        self.logger.info(f"过滤后保留 {len(results)} 条 (score >= {self.min_score})")

        # 清理临时文件
        batch_fasta.unlink()

        if results:
            df = pd.DataFrame(results)
            batch_file = self.output_dir / "batches" / f"batch_{batch_id}_results.tsv"
            df.to_csv(batch_file, sep='\t', index=False)
            return batch_file

        return None

    def run(self):
        """运行完整流程"""
        self.logger.info("\n" + "="*60)
        self.logger.info("Retron RT HMM 搜索流程")
        self.logger.info("="*60)

        start_time = time.time()
        batch_files = []

        for batch in self.batch_plan:
            batch_file = self._process_batch(batch)
            if batch_file:
                batch_files.append(batch_file)
            gc.collect()

        # 合并结果
        self.logger.info("\n合并批次结果...")

        if not batch_files:
            self.logger.warning("未找到任何匹配")
            return False

        dfs = []
        for bf in batch_files:
            dfs.append(pd.read_csv(bf, sep='\t'))

        final_df = pd.concat(dfs, ignore_index=True)

        # 去重：每个位置保留最高分
        before = len(final_df)
        final_df = final_df.sort_values('hmm_score', ascending=False)
        final_df = final_df.drop_duplicates(
            subset=['genome', 'protein_id'],
            keep='first'
        )
        after = len(final_df)

        if before > after:
            self.logger.info(f"去重: {before} -> {after} 条")

        # 保存结果
        results_file = self.output_dir / "reports" / "complete_results.tsv"
        final_df.to_csv(results_file, sep='\t', index=False)
        self.logger.info(f"✓ 完整结果: {results_file}")

        # 生成摘要
        self._generate_summary(final_df, time.time() - start_time)

        # 清理临时目录
        shutil.rmtree(self.temp_dir)

        return True

    def _generate_summary(self, df, elapsed_time):
        """生成摘要报告"""
        summary_file = self.output_dir / "reports" / "summary.txt"

        with open(summary_file, 'w') as f:
            f.write("Retron RT HMM 搜索摘要\n")
            f.write("="*50 + "\n\n")
            f.write(f"HMM模型: {self.hmm_file}\n")
            f.write(f"数据库: {self.database_dir}\n")
            f.write(f"运行时间: {elapsed_time:.1f} 秒\n\n")

            f.write(f"总匹配数: {len(df)}\n")
            f.write(f"涉及基因组: {df['genome'].nunique()}\n\n")

            f.write("HMM Score 分布:\n")
            f.write(f"  最小: {df['hmm_score'].min():.1f}\n")
            f.write(f"  最大: {df['hmm_score'].max():.1f}\n")
            f.write(f"  平均: {df['hmm_score'].mean():.1f}\n")
            f.write(f"  中位数: {df['hmm_score'].median():.1f}\n\n")

            f.write("核心区域长度 (env_to - env_from):\n")
            f.write(f"  最小: {df['env_length'].min()} aa\n")
            f.write(f"  最大: {df['env_length'].max()} aa\n")
            f.write(f"  平均: {df['env_length'].mean():.1f} aa\n")

        self.logger.info(f"✓ 摘要报告: {summary_file}")

        # 打印到终端
        print("\n" + "="*60)
        print("搜索完成!")
        print("="*60)
        print(f"总匹配数: {len(df)}")
        print(f"涉及基因组: {df['genome'].nunique()}")
        print(f"Score范围: {df['hmm_score'].min():.1f} - {df['hmm_score'].max():.1f}")
        print(f"核心区域长度: {df['env_length'].min()}-{df['env_length'].max()} aa")
        print(f"\n输出文件: {self.output_dir / 'reports' / 'complete_results.tsv'}")
        print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description="使用Retron RT HMM模型搜索候选序列",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法
  python step01_hmmsearch.py -m Retron_RT.hmm -d /path/to/faa_dir -o hmm_results

  # 调整参数
  python step01_hmmsearch.py -m Retron_RT.hmm -d /path/to/faa_dir \\
      --min-score 30 --evalue 1e-10 --threads 8

输出文件:
  - reports/complete_results.tsv  : 完整结果 (含env_from/env_to边界)
  - reports/summary.txt           : 摘要统计
        """
    )

    parser.add_argument('-m', '--hmm', required=True,
                        help='HMM模型文件 (Retron_RT.hmm)')
    parser.add_argument('-d', '--database', required=True,
                        help='FAA蛋白序列目录')
    parser.add_argument('-o', '--output', default='hmm_search_results',
                        help='输出目录 (默认: hmm_search_results)')
    parser.add_argument('--batch-size', type=int, default=100,
                        help='批次大小 (默认: 100)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='线程数 (默认: 4)')
    parser.add_argument('-e', '--evalue', type=float, default=1e-5,
                        help='E-value阈值 (默认: 1e-5)')
    parser.add_argument('--min-score', type=float, default=25.0,
                        help='最小domain score (默认: 25.0)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='详细输出')

    args = parser.parse_args()

    pipeline = HMMSearchPipeline(
        hmm_file=args.hmm,
        database_dir=args.database,
        output_dir=args.output,
        batch_size=args.batch_size,
        threads=args.threads,
        evalue=args.evalue,
        min_score=args.min_score
    )

    success = pipeline.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
