#!/usr/bin/env python3
"""
流式处理版本：整合的序列搜索和antiSMASH BGC匹配完整流程
核心改进：每个批次走完整流程（搜索→修复→匹配→保存→释放内存）

优化点：
1. ✅ 每批次完整处理，避免内存累积
2. ✅ sequence_positions和antismash_cache在初始化时加载一次
3. ✅ 批次结果立即处理并保存，不在内存中累积
4. ✅ 最后只需文件级别的合并操作
5. ✅ 大幅降低内存峰值
6. ✅ 使用整个基因位置匹配BGC（更准确）
7. ✅ [新增] 修复Prokka/NCBI ID格式不匹配问题
8. ✅ [新增] 优化比对速度（预计算Map + C语言加速）
9. ✅ [新增] 断点续传功能
"""

import os
import sys
import subprocess
import logging
import shutil
from pathlib import Path
import time
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Align  # [新增] 引入高性能比对模块
import pandas as pd
import argparse
import gc
import psutil
import re
from collections import defaultdict
import tempfile

# 尝试导入intervaltree用于快速区间查询
try:
    from intervaltree import IntervalTree, Interval
    HAS_INTERVALTREE = True
except ImportError:
    HAS_INTERVALTREE = False
    print("⚠️ 建议安装intervaltree以提升性能: pip install intervaltree")

# ============================================================================
# 工具函数部分
# ============================================================================

def get_memory_usage():
    """获取当前内存使用情况（MB）"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

def extract_contig_from_protein_id(protein_id):
    """
    从protein_id中提取contig名称 (增强版，修复Prokka/NCBI格式问题)
    
    示例：
    - NZ_JAIQDG010000100.1_1 → NZ_JAIQDG010000100.1
    - gnl|Prokka|assembly_1_120 → gnl|Prokka|assembly_1
    """
    # 移除可能的 query 描述部分 (空格之后)
    protein_id = protein_id.split()[0]
    
    # [修复] 针对 Prokka 格式: gnl|Prokka|assembly_1_120 -> gnl|Prokka|assembly_1
    if '|' in protein_id:
        parts = protein_id.rsplit('_', 1)
        # 只有当最后一部分是纯数字时，才认为是基因编号并切除
        if len(parts) > 1 and parts[-1].isdigit():
            return parts[0]
            
    # 标准处理
    if '_' in protein_id:
        parts = protein_id.rsplit('_', 1)
        if parts[-1].isdigit():
            return parts[0]
    
    match = re.match(r'([A-Z_]+\d+\.\d+)', protein_id)
    if match:
        return match.group(1)
    
    return protein_id

def detect_sequence_type(file_path):
    """自动检测序列类型：核酸(nucleotide)或蛋白(protein) - 使用流式处理"""
    try:
        sample_count = 0
        total_nucleotide_ratio = 0
        max_samples = 5
        
        for seq in SeqIO.parse(file_path, 'fasta'):
            if sample_count >= max_samples:
                break
            
            seq_str = str(seq.seq[:200]).upper()
            nucleotides = set('ATCGN')
            nucleotide_count = sum(1 for c in seq_str if c in nucleotides)
            ratio = nucleotide_count / len(seq_str) if seq_str else 0
            total_nucleotide_ratio += ratio
            sample_count += 1
        
        if sample_count == 0:
            return "empty"
        
        avg_ratio = total_nucleotide_ratio / sample_count
        seq_type = "nucleotide" if avg_ratio > 0.85 else "protein"
        return seq_type
        
    except Exception as e:
        print(f"警告: 无法检测序列类型 {file_path}: {e}")
        return "unknown"

def determine_search_strategy(query_type, db_type):
    """根据query和database类型自动选择搜索策略"""
    strategies = {
        ('protein', 'protein'): ('blastp', 'blastp', '蛋白-蛋白搜索'),
        ('protein', 'nucleotide'): ('tblastn', None, '蛋白查询-核酸数据库（翻译）'),
        ('nucleotide', 'protein'): ('blastx', 'blastx', '核酸查询（翻译）-蛋白数据库'),
        ('nucleotide', 'nucleotide'): ('blastn', None, '核酸-核酸搜索'),
    }
    
    key = (query_type, db_type)
    if key in strategies:
        return strategies[key]
    else:
        return ('blastp', 'blastp', '默认蛋白搜索')

# ============================================================================
# 位置信息提取和转换
# ============================================================================

def parse_sequence_location(description):
    """从FASTA描述中解析位置信息（核苷酸坐标）"""
    patterns = [
        r'Location[:\s]+([\d,]+)\s*-\s*([\d,]+)',
        r'location[=:]\s*(\d+)\.\.(\d+)',
        r'coordinates[=:]\s*(\d+)-(\d+)',
        r'\[location[=:]\s*(\d+)\.\.(\d+)\]',
        r'#\s*(\d+)\s*#\s*(\d+)\s*#',
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
    
    if re.search(r'strand[=:]?\s*-', description, re.IGNORECASE):
        return '-'
    elif re.search(r'strand[=:]?\s*\+', description, re.IGNORECASE):
        return '+'
    
    return None

def extract_positions_from_fasta(fasta_file):
    """从FAA文件中提取所有序列的位置信息"""
    positions = {}

    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq_id = record.id
        description = record.description
        contig_name = extract_contig_from_protein_id(seq_id)
        
        nt_start, nt_end = parse_sequence_location(description)
        strand = parse_strand_from_description(description)
        
        if strand is None and nt_start and nt_end:
            strand = '+' if nt_start < nt_end else '-'
        
        positions[seq_id] = {
            'nt_start': nt_start,
            'nt_end': nt_end,
            'strand': strand,
            'aa_length': len(record.seq),
            'description': description,
            'contig': contig_name,
        }
    
    return positions

def scan_database_directory_for_positions(database_dir, logger):
    """扫描数据库目录，提取所有序列的位置信息"""
    database_dir = Path(database_dir)
    all_positions = {}
    
    fasta_extensions = ['*.faa', '*.fasta', '*.fa', '*.pep']
    fasta_files = []
    
    for ext in fasta_extensions:
        fasta_files.extend(database_dir.rglob(ext))
    
    logger.info(f"扫描位置信息: 找到 {len(fasta_files)} 个序列文件")
    
    for idx, fasta_file in enumerate(fasta_files):
        genome_name = fasta_file.stem
        
        try:
            positions = extract_positions_from_fasta(fasta_file)
            with_location = sum(1 for p in positions.values() if p['nt_start'])
            
            if with_location > 0:
                logger.info(f"  ✓ {genome_name}: {with_location}/{len(positions)} 个序列有位置信息")
            
            for seq_id, pos_info in positions.items():
                full_id = f"{genome_name}___{seq_id}"
                all_positions[full_id] = pos_info
                
                if seq_id not in all_positions:
                    all_positions[seq_id] = pos_info
            
            if (idx + 1) % 10 == 0:
                gc.collect()
                mem_mb = get_memory_usage()
                logger.info(f"  💾 已处理 {idx + 1}/{len(fasta_files)} 个文件，当前内存: {mem_mb:.1f} MB")
                
        except Exception as e:
            logger.warning(f"  ✗ 处理 {genome_name} 失败: {e}")
    
    logger.info(f"成功提取 {len(all_positions)} 个序列的位置信息")
    return all_positions

def aa_to_nt_position(aa_pos, nt_start, nt_end, strand):
    """将氨基酸位置转换为核苷酸位置（保留此函数以防其他地方调用，但不再使用）"""
    if strand == '+':
        return nt_start + (aa_pos - 1) * 3
    else:
        return nt_end - (aa_pos - 1) * 3

# ============================================================================
# antiSMASH结果解析 (Global Class, mostly used for reference structure)
# ============================================================================

class AntiSMASHCache:
    """antiSMASH数据缓存类，使用interval tree加速查询"""
    
    def __init__(self, antismash_dir, logger):
        self.logger = logger
        self.data = defaultdict(lambda: defaultdict(lambda: {'regions': [], 'cds': []}))
        self.interval_trees = {}
        self._load_data(antismash_dir)
        self._build_interval_trees()
    
    def _load_data(self, antismash_dir):
        """加载antiSMASH数据"""
        antismash_dir = Path(antismash_dir)
        genome_dirs = [d for d in antismash_dir.iterdir() if d.is_dir()]
        if not genome_dirs:
            genome_dirs = [antismash_dir]
        
        self.logger.info(f"   扫描 {len(genome_dirs)} 个antiSMASH目录...")
        
        total_regions = 0
        total_cds = 0
        
        for genome_dir in genome_dirs:
            genome_name = genome_dir.name
            gbk_files = [f for f in genome_dir.glob("*.gbk") if 'region' not in f.name.lower()]
            
            if not gbk_files:
                continue
            
            for gbk_file in gbk_files:
                try:
                    for record in SeqIO.parse(gbk_file, 'genbank'):
                        contig_id = record.id
                        
                        for feature in record.features:
                            if feature.type == 'region' and 'region_number' in feature.qualifiers:
                                start = int(feature.location.start) + 1
                                end = int(feature.location.end)
                                region_num = feature.qualifiers['region_number'][0]
                                
                                region_type = 'Unknown'
                                if 'product' in feature.qualifiers:
                                    products = feature.qualifiers['product']
                                    region_type = ';'.join(products) if isinstance(products, list) else products
                                
                                self.data[genome_name][contig_id]['regions'].append({
                                    'start': start,
                                    'end': end,
                                    'type': region_type,
                                    'region_number': region_num
                                })
                                total_regions += 1
                        
                        for feature in record.features:
                            if feature.type == 'CDS':
                                start = int(feature.location.start) + 1
                                end = int(feature.location.end)
                                strand = '+' if feature.location.strand == 1 else '-'
                                
                                locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                                protein_id = feature.qualifiers.get('protein_id', [''])[0]
                                product = feature.qualifiers.get('product', [''])[0]
                                gene_kind = feature.qualifiers.get('gene_kind', [''])[0]
                                
                                self.data[genome_name][contig_id]['cds'].append({
                                    'start': start,
                                    'end': end,
                                    'strand': strand,
                                    'locus_tag': locus_tag,
                                    'protein_id': protein_id,
                                    'product': product,
                                    'gene_kind': gene_kind
                                })
                                total_cds += 1
                
                except Exception as e:
                    self.logger.warning(f"   解析 {gbk_file.name} 失败: {e}")
        
        self.logger.info(f"   提取了 {total_regions} 个regions, {total_cds} 个CDS")
    
    def _build_interval_trees(self):
        """构建interval tree用于快速区间查询"""
        self.logger.info("   ⚠️ 跳过interval tree构建，使用线性搜索")
        return
    
    def find_overlapping_cds(self, genome, contig, start, end):
        """使用interval tree快速查找重叠的CDS"""
        if HAS_INTERVALTREE and genome in self.interval_trees:
            if contig in self.interval_trees[genome]:
                tree = self.interval_trees[genome][contig]
                overlaps = tree.overlap(start, end + 1)
                return [interval.data for interval in overlaps]
        
        if genome in self.data and contig in self.data[genome]:
            cds_list = self.data[genome][contig]['cds']
            result = []
            for cds in cds_list:
                if not (end < cds['start'] or start > cds['end']):
                    result.append(cds)
            return result
        
        return []

# ============================================================================
# BLAST/DIAMOND搜索函数
# ============================================================================

def parse_blast_hit(title, genome_mapping):
    """解析BLAST结果中的基因组和序列信息"""
    if title in genome_mapping:
        return genome_mapping[title]['genome'], genome_mapping[title]['sequence_id']
    
    for mapped_id, info in genome_mapping.items():
        if title in mapped_id or mapped_id in title:
            return info['genome'], info['sequence_id']
    
    if "___" in title:
        parts = title.split("___")
        if len(parts) >= 2:
            return parts[0].strip(), parts[1].split()[0]
    
    return title.strip(), title.strip()

def run_diamond_search(query_name, query_file, db_path, output_dir, genome_mapping, 
                      diamond_mode='blastp', threads=4):
    """运行DIAMOND搜索"""
    output_file = output_dir / "diamond" / f"{query_name}_diamond.tsv"
    
    if not shutil.which('diamond'):
        return []
    
    cmd = [
        'diamond', diamond_mode,
        '--query', str(query_file),
        '--db', db_path,
        '--out', str(output_file),
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 
                   'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen',
        '--evalue', '1e-5',
        '--max-target-seqs', '500',
        '--threads', str(threads),
        '--sensitive'
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=1800)
        
        hits = []
        with open(output_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 12:
                    try:
                        query_id = fields[0]
                        subject_id = fields[1]
                        identity = float(fields[2])
                        align_length = int(fields[3])
                        evalue = float(fields[4])
                        bitscore = float(fields[5])
                        qstart = int(fields[6])
                        qend = int(fields[7])
                        sstart = int(fields[8])
                        send = int(fields[9])
                        qlen = int(fields[10])
                        slen = int(fields[11])
                        
                        coverage = (align_length / qlen) * 100
                        strand = '+' if sstart < send else '-'
                        
                        genome_name, seq_id = parse_blast_hit(subject_id, genome_mapping)
                        
                        hits.append({
                            'query': query_name,
                            'query_id': query_id,
                            'genome': genome_name,
                            'contig': seq_id,
                            'position': f"{min(sstart, send)}-{max(sstart, send)}",
                            'start': min(sstart, send),
                            'end': max(sstart, send),
                            'strand': strand,
                            'identity': identity,
                            'coverage': coverage,
                            'evalue': evalue,
                            'bitscore': bitscore,
                            'alignment_length': align_length,
                            'query_start': qstart,
                            'query_end': qend,
                            'subject_start': sstart,
                            'subject_end': send,
                            'query_length': qlen,
                            'subject_length': slen,
                            'method': f'DIAMOND_{diamond_mode.upper()}'
                        })
                    except (ValueError, IndexError):
                        continue
        
        return hits
    except Exception as e:
        return []

def run_blast_search(query_name, query_file, db_path, output_dir, genome_mapping, 
                     blast_mode='blastp', threads=2):
    """运行BLAST搜索"""
    output_file = output_dir / "blast" / f"{query_name}_{blast_mode}.xml"
    
    blast_commands = {
        'blastn': 'blastn',
        'blastp': 'blastp',
        'tblastn': 'tblastn',
        'blastx': 'blastx'
    }
    
    cmd = [
        blast_commands[blast_mode],
        '-query', str(query_file),
        '-db', db_path,
        '-out', str(output_file),
        '-outfmt', '5',
        '-evalue', '1e-5',
        '-max_target_seqs', '500',
        '-num_threads', str(threads)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=3600)
        
        hits = []
        with open(output_file) as f:
            for blast_record in NCBIXML.parse(f):
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        identity = (hsp.identities / hsp.align_length) * 100
                        coverage = (hsp.align_length / blast_record.query_length) * 100
                        
                        genome_name, seq_id = parse_blast_hit(alignment.title, genome_mapping)
                        strand = '+' if hsp.sbjct_start < hsp.sbjct_end else '-'
                        
                        hits.append({
                            'query': query_name,
                            'query_id': blast_record.query,
                            'genome': genome_name,
                            'contig': seq_id,
                            'position': f"{min(hsp.sbjct_start, hsp.sbjct_end)}-{max(hsp.sbjct_start, hsp.sbjct_end)}",
                            'start': min(hsp.sbjct_start, hsp.sbjct_end),
                            'end': max(hsp.sbjct_start, hsp.sbjct_end),
                            'strand': strand,
                            'identity': identity,
                            'coverage': coverage,
                            'evalue': hsp.expect,
                            'bitscore': hsp.bits,
                            'alignment_length': hsp.align_length,
                            'query_start': hsp.query_start,
                            'query_end': hsp.query_end,
                            'subject_start': hsp.sbjct_start,
                            'subject_end': hsp.sbjct_end,
                            'query_length': blast_record.query_length,
                            'subject_length': alignment.length,
                            'method': f'BLAST_{blast_mode.upper()}'
                        })
        
        return hits
    except Exception as e:
        return []

# ============================================================================
# 🆕 流式处理主类 - 每批次完整处理
# ============================================================================

class StreamProcessingPipeline:
    """流式处理版本：每批次走完整流程"""
    
    def __init__(self, database_dir, query_dir, output_dir="search_results",
                 batch_size=100, threads=4, use_diamond=True, use_blast=False,
                 min_identity=30, min_coverage=50, antismash_dir=None):
        
        self.database_dir = Path(database_dir)
        self.query_dir = Path(query_dir)
        self.output_dir = Path(output_dir)
        self.antismash_dir = Path(antismash_dir) if antismash_dir else None
        self.batch_size = batch_size
        self.threads = threads
        self.use_diamond = use_diamond
        self.use_blast = use_blast
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        
        self.temp_dir = None
        
        self._create_output_structure()
        self._setup_logging()
        self._check_dependencies()
        
        self.logger.info("📋 步骤0: 分析数据类型和预加载...")
        self.database_info = self._analyze_database_files()
        self.query_info = self._analyze_query_files()
        self._determine_search_strategy()
        
        # 🆕 关键改进：不在初始化时预加载，而是每批次加载对应的数据
        self.logger.info("\n" + "="*80)
        self.logger.info("📦 批次加载模式（每批次只加载需要的数据）")
        self.logger.info("="*80)
        
        # 不预加载，每批次动态加载
        self.sequence_positions = None
        self.antismash_cache = None
        
        self.logger.info(f"💾 初始化后内存: {get_memory_usage():.1f} MB")
        
        self.batch_plan = self._create_batch_plan()
    
    def _create_output_structure(self):
        """创建输出目录结构"""
        directories = ["diamond", "blast", "reports", "logs", "batches"]
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for dir_name in directories:
            (self.output_dir / dir_name).mkdir(parents=True, exist_ok=True)
        
        self.temp_dir = Path(tempfile.mkdtemp(prefix="blast_search_"))
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"临时文件目录: {self.temp_dir}")
    
    def _setup_logging(self):
        """设置日志系统"""
        log_file = self.output_dir / "logs" / "pipeline.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def _check_dependencies(self):
        """检查必要的依赖工具"""
        self.diamond_available = shutil.which('diamond') is not None
        self.blast_available = shutil.which('blastp') is not None
        
        if self.use_diamond and not self.diamond_available:
            self.logger.warning("DIAMOND不可用，将禁用DIAMOND功能")
            self.use_diamond = False
        
        if self.use_blast and not self.blast_available:
            self.logger.warning("BLAST不可用，将禁用BLAST功能")
            self.use_blast = False
        
        if not self.use_diamond and not self.use_blast:
            self.logger.error("DIAMOND和BLAST都不可用，无法进行搜索")
            sys.exit(1)
    
    def _analyze_database_files(self):
        """分析数据库文件并自动检测类型"""
        self.logger.info("📚 分析数据库文件...")
        
        db_files = list(self.database_dir.glob("*.fasta")) + \
                   list(self.database_dir.glob("*.fa")) + \
                   list(self.database_dir.glob("*.fna")) + \
                   list(self.database_dir.glob("*.faa")) + \
                   list(self.database_dir.glob("*.pep"))
        
        if not db_files:
            db_files = list(self.database_dir.rglob("*.fasta")) + \
                       list(self.database_dir.rglob("*.fa")) + \
                       list(self.database_dir.rglob("*.fna")) + \
                       list(self.database_dir.rglob("*.faa")) + \
                       list(self.database_dir.rglob("*.pep"))
        
        if not db_files:
            self.logger.error(f"未在 {self.database_dir} 中找到序列文件")
            sys.exit(1)
        
        database_info = {}
        db_types = []
        
        for db_file in db_files:
            base_name = db_file.stem
            seq_type = detect_sequence_type(db_file)
            
            database_info[base_name] = {
                'file': db_file,
                'type': seq_type
            }
            
            if seq_type != "unknown":
                db_types.append(seq_type)
        
        if db_types:
            db_type_counts = pd.Series(db_types).value_counts()
            self.database_type = db_type_counts.index[0]
            self.logger.info(f"📊 数据库类型: {self.database_type.upper()}")
        else:
            self.logger.error("无法确定数据库类型")
            sys.exit(1)
        
        return database_info
    
    def _analyze_query_files(self):
        """分析查询序列文件并自动检测类型"""
        self.logger.info("🎯 分析查询序列...")
        
        query_files = list(self.query_dir.glob("*.fasta")) + \
                     list(self.query_dir.glob("*.fa")) + \
                     list(self.query_dir.glob("*.fna")) + \
                     list(self.query_dir.glob("*.faa")) + \
                     list(self.query_dir.glob("*.pep"))
        
        if not query_files:
            self.logger.error(f"未在 {self.query_dir} 中找到查询序列文件")
            sys.exit(1)
        
        query_info = {}
        query_types = []
        
        for query_file in query_files:
            query_name = query_file.stem
            seq_type = detect_sequence_type(query_file)
            seq_count = sum(1 for _ in SeqIO.parse(query_file, 'fasta'))
            
            query_info[query_name] = {
                'file': query_file,
                'type': seq_type,
                'sequence_count': seq_count
            }
            
            if seq_type != "unknown":
                query_types.append(seq_type)
        
        if query_types:
            query_type_counts = pd.Series(query_types).value_counts()
            self.query_type = query_type_counts.index[0]
            self.logger.info(f"📊 查询类型: {self.query_type.upper()}")
        else:
            self.logger.error("无法确定查询序列类型")
            sys.exit(1)
        
        return query_info
    
    def _determine_search_strategy(self):
        """根据数据类型自动确定搜索策略"""
        self.logger.info("🤖 确定搜索策略...")
        
        blast_mode, diamond_mode, description = determine_search_strategy(
            self.query_type, self.database_type
        )
        
        self.blast_mode = blast_mode
        self.diamond_mode = diamond_mode
        self.is_protein_search = (self.query_type == 'protein')
        
        self.logger.info(f"   策略: {description}")
        self.logger.info(f"   BLAST模式: {blast_mode.upper()}")
        
        if diamond_mode:
            self.logger.info(f"   DIAMOND模式: {diamond_mode.upper()}")
        else:
            if self.use_diamond and not self.use_blast:
                self.logger.warning("   切换到BLAST模式")
                self.use_blast = True
                self.use_diamond = False
    
    def _create_batch_plan(self):
        """创建批次处理计划"""
        db_names = list(self.database_info.keys())
        
        batches = []
        for i in range(0, len(db_names), self.batch_size):
            batch_items = db_names[i:i + self.batch_size]
            batches.append({
                'batch_id': i // self.batch_size + 1,
                'items': batch_items,
                'size': len(batch_items)
            })
        
        self.logger.info(f"📋 创建 {len(batches)} 个批次 (每批 ~{self.batch_size} 个文件)")
        
        return batches
    
    def _create_batch_database(self, batch_items, batch_id):
        """创建批次数据库"""
        batch_seq_file = self.temp_dir / f"batch_{batch_id}_db.fasta"
        genome_mapping = {}
        
        seq_count = 0
        with open(batch_seq_file, 'w') as out_f:
            for item_name in batch_items:
                info = self.database_info[item_name]
                seq_file = info['file']
                
                for record in SeqIO.parse(seq_file, 'fasta'):
                    seq_count += 1
                    original_id = record.id
                    new_id = f"{item_name}___{original_id}"
                    
                    record.id = new_id
                    record.description = f"{item_name} {record.description}"
                    SeqIO.write(record, out_f, 'fasta')
                    
                    genome_mapping[new_id] = {
                        'genome': item_name,
                        'sequence_id': original_id
                    }
        
        db_path = None
        
        if self.use_diamond and self.diamond_mode:
            diamond_db_path = self.temp_dir / f"batch_{batch_id}_diamond_db"
            cmd = ['diamond', 'makedb', '--in', str(batch_seq_file), '--db', str(diamond_db_path)]
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                db_path = str(diamond_db_path)
            except Exception as e:
                self.logger.warning(f"   DIAMOND数据库创建失败: {e}")
        
        if self.use_blast or not db_path:
            blast_db_path = self.temp_dir / f"batch_{batch_id}_blast_db"
            dbtype = 'nucl' if self.database_type == 'nucleotide' else 'prot'
            cmd = ['makeblastdb', '-in', str(batch_seq_file), '-dbtype', dbtype, 
                   '-out', str(blast_db_path)]
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                if not db_path:
                    db_path = str(blast_db_path)
            except Exception as e:
                self.logger.error(f"   BLAST数据库创建失败: {e}")
                if not db_path:
                    raise
        
        return db_path, genome_mapping, batch_seq_file
    
    def _search_batch(self, db_path, genome_mapping, batch_id):
        """搜索单个批次"""
        batch_results = []
        
        for query_name, query_info in self.query_info.items():
            query_file = query_info['file']
            
            if self.use_diamond and self.diamond_mode:
                hits = run_diamond_search(
                    query_name, query_file, str(db_path),
                    self.output_dir, genome_mapping, 
                    self.diamond_mode, self.threads
                )
                batch_results.extend(hits)
                self.logger.info(f"   DIAMOND {query_name}: {len(hits)} 个匹配")
            
            if self.use_blast:
                hits = run_blast_search(
                    query_name, query_file, str(db_path),
                    self.output_dir, genome_mapping, 
                    self.blast_mode, self.threads
                )
                batch_results.extend(hits)
                self.logger.info(f"   BLAST {query_name}: {len(hits)} 个匹配")
        
        return batch_results
    
    def _cleanup_batch_files(self, batch_seq_file, db_path, batch_id):
        """清理批次临时文件 (含内存回收优化)"""
        try:
            if batch_seq_file.exists():
                batch_seq_file.unlink()
            
            for ext in ['', '.dmnd', '.pin', '.phr', '.psq', '.nhr', '.nin', '.nsq']:
                db_file = Path(str(db_path) + ext)
                if db_file.exists():
                    db_file.unlink()
            
            # [优化] 显式垃圾回收，处理 Biopython 可能残留的对象引用
            gc.collect()
            
            self.logger.info(f"   ✓ 清理批次 {batch_id} 临时文件与内存")
        except Exception as e:
            self.logger.warning(f"   清理临时文件失败: {e}")
    
    def _load_batch_positions(self, batch_items):
        """加载批次涉及的基因组位置信息"""
        if not self.is_protein_search:
            return {}
        
        self.logger.info("   📍 加载批次位置信息...")
        batch_positions = {}
        
        for item_name in batch_items:
            info = self.database_info[item_name]
            fasta_file = info['file']
            
            try:
                positions = extract_positions_from_fasta(fasta_file)
                
                for seq_id, pos_info in positions.items():
                    full_id = f"{item_name}___{seq_id}"
                    batch_positions[full_id] = pos_info
                    
                    if seq_id not in batch_positions:
                        batch_positions[seq_id] = pos_info
                        
            except Exception as e:
                self.logger.warning(f"   提取 {item_name} 位置信息失败: {e}")
        
        self.logger.info(f"   ✓ 加载了 {len(batch_positions)} 条位置信息")
        return batch_positions
    
    def _load_batch_antismash(self, batch_items):
        """加载批次涉及的基因组antiSMASH数据"""
        # 调试信息
        self.logger.info(f"   📍 检查antiSMASH: dir={self.antismash_dir}, exists={self.antismash_dir.exists() if self.antismash_dir else 'None'}")
        
        if not self.antismash_dir or not self.antismash_dir.exists():
            self.logger.info("   ⚠️ antiSMASH目录不可用，跳过BGC匹配")
            return None
        
        self.logger.info("   🧬 加载批次antiSMASH数据...")
        
        # 创建临时的batch级别的antismash cache
        class BatchAntiSMASHCache:
            def __init__(self, antismash_dir, batch_items, logger):
                self.logger = logger
                self.data = defaultdict(lambda: defaultdict(lambda: {'regions': [], 'cds': []}))
                # [优化] 新增：预计算 Contig 名称映射表 {genome: {normalized_name: real_name}}
                self.contig_map = defaultdict(dict)
                self.interval_trees = {}
                self._load_batch_data(antismash_dir, batch_items)
                self._build_interval_trees()
            
            def _load_batch_data(self, antismash_dir, batch_items):
                """只加载当前批次基因组的antiSMASH数据"""
                total_regions = 0
                total_cds = 0
                found_genomes = 0
                
                self.logger.info(f"   查找antiSMASH数据：批次包含 {len(batch_items)} 个基因组")

                # Bug修复：预索引antiSMASH目录，支持核心数字匹配
                def extract_core_number(name):
                    """提取核心数字用于模糊匹配"""
                    match = re.search(r'(\d{9,}(?:\.\d+)?)', str(name))
                    return match.group(1) if match else str(name)

                antismash_dirs_map = {}
                for d in antismash_dir.iterdir():
                    if d.is_dir():
                        antismash_dirs_map[d.name] = d
                        antismash_dirs_map[extract_core_number(d.name)] = d

                for genome_name in batch_items:
                    # 尝试直接匹配
                    genome_dir = antismash_dir / genome_name

                    # 如果直接匹配失败，尝试核心数字匹配
                    if not genome_dir.exists():
                        core_num = extract_core_number(genome_name)
                        if core_num in antismash_dirs_map:
                            genome_dir = antismash_dirs_map[core_num]
                            self.logger.debug(f"   ✓ 核心数字匹配: {genome_name} -> {genome_dir.name}")

                    if not genome_dir.exists():
                        self.logger.debug(f"   ⚠️ 未找到目录: {genome_dir}")
                        continue

                    if not genome_dir.is_dir():
                        self.logger.debug(f"   ⚠️ 不是目录: {genome_dir}")
                        continue
                    
                    found_genomes += 1
                    
                    gbk_files = [f for f in genome_dir.glob("*.gbk") if 'region' not in f.name.lower()]
                    
                    if not gbk_files:
                        self.logger.debug(f"   ⚠️ {genome_name}: 未找到GBK文件")
                        continue
                    
                    self.logger.debug(f"   ✓ {genome_name}: 找到 {len(gbk_files)} 个GBK文件")
                    
                    for gbk_file in gbk_files:
                        try:
                            # 添加文件检查
                            if not gbk_file.exists() or gbk_file.stat().st_size == 0:
                                self.logger.warning(f"   ⚠️ {gbk_file.name}: 文件不存在或为空")
                                continue

                            # 直接处理所有记录，不限制数量
                            try:
                                for record in SeqIO.parse(gbk_file, 'genbank'):

                                    try:
                                        contig_id = record.id

                                        # [优化] 预先计算并存储 Contig 名称映射
                                        # 将 "gnlProkkaassembly_1" -> "gnlprokkaassembly1"
                                        norm_name = re.sub(r'[^a-zA-Z0-9]', '', str(contig_id)).lower()
                                        self.contig_map[genome_name][norm_name] = contig_id

                                        # 处理regions
                                        for feature in record.features:
                                            try:
                                                if feature.type == 'region' and 'region_number' in feature.qualifiers:
                                                    start = int(feature.location.start) + 1
                                                    end = int(feature.location.end)
                                                    region_num = feature.qualifiers['region_number'][0]

                                                    region_type = 'Unknown'
                                                    if 'product' in feature.qualifiers:
                                                        products = feature.qualifiers['product']
                                                        region_type = ';'.join(products) if isinstance(products, list) else products

                                                    self.data[genome_name][contig_id]['regions'].append({
                                                        'start': start,
                                                        'end': end,
                                                        'type': region_type,
                                                        'region_number': region_num
                                                    })
                                                    total_regions += 1
                                            except Exception as e:
                                                # 跳过有问题的feature，继续处理下一个
                                                continue

                                        # 处理CDS
                                        for feature in record.features:
                                            try:
                                                if feature.type == 'CDS':
                                                    start = int(feature.location.start) + 1
                                                    end = int(feature.location.end)
                                                    strand = '+' if feature.location.strand == 1 else '-'

                                                    locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                                                    protein_id = feature.qualifiers.get('protein_id', [''])[0]
                                                    product = feature.qualifiers.get('product', [''])[0]
                                                    gene_kind = feature.qualifiers.get('gene_kind', [''])[0]

                                                    self.data[genome_name][contig_id]['cds'].append({
                                                        'start': start,
                                                        'end': end,
                                                        'strand': strand,
                                                        'locus_tag': locus_tag,
                                                        'protein_id': protein_id,
                                                        'product': product,
                                                        'gene_kind': gene_kind
                                                    })
                                                    total_cds += 1
                                            except Exception as e:
                                                # 跳过有问题的CDS feature，继续处理下一个
                                                continue

                                    except Exception as e:
                                        # 跳过有问题的record，继续处理下一个
                                        continue

                            except MemoryError:
                                self.logger.error(f"   ❌ {gbk_file.name}: 内存不足，跳过该文件")
                                gc.collect()
                                continue
                            except KeyboardInterrupt:
                                self.logger.error(f"   ⚠️ {gbk_file.name}: 用户中断处理")
                                raise
                            except Exception as e:
                                # 捕获SeqIO解析时的其他错误（包括段错误）
                                self.logger.error(f"   ❌ {gbk_file.name}: SeqIO解析失败 - {e}")
                                # 如果是段错误，尝试清理内存
                                if "segmentation" in str(e).lower() or "segfault" in str(e).lower():
                                    self.logger.warning(f"   ⚠️ 检测到段错误，清理内存并跳过文件 {gbk_file.name}")
                                    gc.collect()
                                continue

                        except Exception as e:
                            self.logger.warning(f"   ❌ 处理 {gbk_file.name} 时发生未预期错误: {e}")
                            gc.collect()
                            continue
                
                self.logger.info(f"   ✓ 匹配到 {found_genomes}/{len(batch_items)} 个基因组目录")
                self.logger.info(f"   ✓ 提取了 {total_regions} 个regions, {total_cds} 个CDS")
            
            def _build_interval_trees(self):
                """构建interval tree用于快速区间查询"""
                # 禁用interval tree以避免内存问题，改用线性搜索
                # 当数据量过大时，interval tree会消耗大量内存导致段错误
                self.logger.info("   ⚠️ 跳过interval tree构建，使用线性搜索")
                return

                # 原代码（已注释）：
                # if not HAS_INTERVALTREE:
                #     return
                #
                # for genome_name, contigs in self.data.items():
                #     self.interval_trees[genome_name] = {}
                #     for contig_name, data in contigs.items():
                #         tree = IntervalTree()
                #
                #         for cds in data['cds']:
                #             tree.addi(cds['start'], cds['end'] + 1, cds)
                #
                #         self.interval_trees[genome_name][contig_name] = tree
            
            def find_overlapping_cds(self, genome, contig, start, end):
                """使用interval tree快速查找重叠的CDS"""
                if HAS_INTERVALTREE and genome in self.interval_trees:
                    if contig in self.interval_trees[genome]:
                        tree = self.interval_trees[genome][contig]
                        overlaps = tree.overlap(start, end + 1)
                        return [interval.data for interval in overlaps]
                
                if genome in self.data and contig in self.data[genome]:
                    cds_list = self.data[genome][contig]['cds']
                    result = []
                    for cds in cds_list:
                        if not (end < cds['start'] or start > cds['end']):
                            result.append(cds)
                    return result
                
                return []
        
        batch_cache = BatchAntiSMASHCache(self.antismash_dir, batch_items, self.logger)
        return batch_cache
    
    def _fix_positions_batch(self, batch_df, batch_positions):
        """修复批次结果的位置信息（使用整个基因的位置，不进行氨基酸到核苷酸的转换）"""
        if not self.is_protein_search or batch_df.empty:
            return batch_df
        
        self.logger.info("   🔧 修复位置...")
        
        batch_df['genome_nt_start'] = None
        batch_df['genome_nt_end'] = None
        batch_df['gene_nt_location'] = None
        batch_df['gene_strand'] = None
        batch_df['position_fixed'] = False
        
        fixed_count = 0
        
        for row in batch_df.itertuples():
            idx = row.Index
            contig = row.contig
            genome = row.genome
            
            possible_ids = [
                f"{genome}___{contig}",
                contig,
            ]
            
            pos_info = None
            for test_id in possible_ids:
                if test_id in batch_positions:
                    pos_info = batch_positions[test_id]
                    break
            
            if pos_info and pos_info['nt_start']:
                # ✅ 直接使用整个基因的位置，不需要转换！
                nt_start = pos_info['nt_start']
                nt_end = pos_info['nt_end']
                strand = pos_info['strand'] or '+'
                
                # 确保起始<结束（升序排列）
                genome_nt_start = min(nt_start, nt_end)
                genome_nt_end = max(nt_start, nt_end)
                
                batch_df.at[idx, 'genome_nt_start'] = genome_nt_start
                batch_df.at[idx, 'genome_nt_end'] = genome_nt_end
                batch_df.at[idx, 'gene_nt_location'] = f"{genome_nt_start}-{genome_nt_end}"
                batch_df.at[idx, 'gene_strand'] = strand
                batch_df.at[idx, 'position_fixed'] = True
                
                fixed_count += 1
        
        self.logger.info(f"   ✓ 修复 {fixed_count}/{len(batch_df)} 条")
        
        return batch_df
    
    def _extract_id_from_fasta(self, record_description):
        """从FAA记录描述中提取ID字段（如ID=1_1）"""
        # 查找 ID=xxx 格式
        id_match = re.search(r'ID=([^|\s]+)', record_description)
        if id_match:
            return id_match.group(1).strip()

        # 查找 GeneID=xxx 格式
        geneid_match = re.search(r'GeneID=([^|\s]+)', record_description)
        if geneid_match:
            return geneid_match.group(1).strip()

        return None

    def _convert_id_to_locus_tag(self, id_str, contig_name):
        """将ID字段转换为locus_tag格式（增强版）"""
        if not id_str:
            return None

        # 情况1：ID=1_1 格式
        if '_' in id_str and id_str.replace('_', '').isdigit():
            # 获取 Contig 基础名
            contig_base = extract_contig_from_protein_id(contig_name)
            
            # 如果 contig_base 包含 assembly_1 这种复杂名称，直接拼接可能会错
            # 但为了保持兼容性，先尝试返回标准拼接
            # 后续的 fuzzy match 会处理具体差异
            return f"{contig_base}_{id_str}"

        # 情况2：已经是完整格式，直接返回
        if id_str.startswith('ctg') or id_str.startswith('contig') or 'Prokka' in id_str:
            return id_str

        # 情况3：尝试从其他格式提取
        if '_' in id_str:
            parts = id_str.split('_')
            if len(parts) >= 2:
                contig_base = extract_contig_from_protein_id(contig_name)
                return f"{contig_base}_{parts[-1]}"

        return id_str

    def _load_gbk_sequences(self, genome, contig_name, force_all_contigs=False):
        """加载GBK文件的序列信息用于序列比对

        Args:
            genome: 基因组名称
            contig_name: contig名称
            force_all_contigs: 是否强制加载所有contig（用于contig匹配失败时）
        """
        gbk_sequences = {}

        if genome not in self.antismash_cache.data:
            return gbk_sequences

        # 加载该genome下的所有contig（当contig_name为None或force_all_contigs=True时）
        contigs_to_search = []
        if force_all_contigs or not contig_name or contig_name not in self.antismash_cache.data[genome]:
            contigs_to_search = list(self.antismash_cache.data[genome].keys())
            self.logger.debug(f"   🔍 全库序列比对: 在 {len(contigs_to_search)} 个contig中搜索")
        else:
            contigs_to_search = [contig_name]

        # 尝试找到对应的GBK文件
        if self.antismash_dir:
            genome_dir = self.antismash_dir / genome
            if genome_dir.exists():
                gbk_files = [f for f in genome_dir.glob("*.gbk") if 'region' not in f.name.lower()]

                for gbk_file in gbk_files:
                    try:
                        # 添加文件检查
                        if not gbk_file.exists() or gbk_file.stat().st_size == 0:
                            continue

                        try:
                            for record in SeqIO.parse(gbk_file, 'genbank'):

                                try:
                                    # 检查contig是否在搜索列表中
                                    if record.id in contigs_to_search or force_all_contigs:
                                        # 提取contig信息
                                        record_contig = record.id

                                        for feature in record.features:
                                            try:
                                                if feature.type == 'CDS':
                                                    locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                                                    if locus_tag:
                                                        # 提取蛋白质序列
                                                        translation = feature.qualifiers.get('translation', [''])[0]
                                                        if translation:
                                                            gbk_sequences[locus_tag] = {
                                                                'sequence': str(translation),
                                                                'start': int(feature.location.start) + 1,
                                                                'end': int(feature.location.end),
                                                                'strand': '+' if feature.location.strand == 1 else '-',
                                                                'contig': record_contig,
                                                                'locus_tag': locus_tag,  # Bug修复：添加locus_tag键
                                                                'product': feature.qualifiers.get('product', [''])[0]
                                                            }
                                            except Exception as e:
                                                # 跳过有问题的feature
                                                continue

                                except Exception as e:
                                    # 跳过有问题的record
                                    continue

                        except MemoryError:
                            self.logger.debug(f"   内存不足，跳过文件: {gbk_file.name}")
                            gc.collect()
                            continue
                        except Exception as e:
                            # 捕获SeqIO解析错误，包括段错误
                            self.logger.debug(f"   SeqIO解析失败: {gbk_file.name}, {e}")
                            if "segmentation" in str(e).lower() or "segfault" in str(e).lower():
                                self.logger.warning(f"   检测到段错误，跳过文件: {gbk_file.name}")
                                gc.collect()
                            continue

                    except Exception as e:
                        self.logger.debug(f"   解析GBK文件失败: {gbk_file.name}, {e}")
                        gc.collect()
                        continue

        return gbk_sequences

    def _get_query_sequence(self, query_id):
        """根据query_id获取查询序列"""
        try:
            # 从query_info中查找对应的文件
            for query_name, query_info in self.query_info.items():
                query_file = query_info['file']
                for record in SeqIO.parse(query_file, 'fasta'):
                    if record.id == query_id:
                        return str(record.seq)
            return None
        except Exception as e:
            self.logger.debug(f"   获取查询序列失败 {query_id}: {e}")
            return None

    def _perform_sequence_alignment(self, query_seq, gbk_sequences):
        """执行序列比对，返回最佳匹配

        Args:
            query_seq: 查询序列
            gbk_sequences: GBK中的序列字典

        Returns:
            dict: 包含最佳匹配的字典，格式为 {
                'cds': cds_info,
                'locus_tag': locus_tag,
                'similarity': similarity_score,
                'contig': contig_name
            }
        """
        best_match = None
        best_score = 0

        for locus_tag, seq_info in gbk_sequences.items():
            target_seq = seq_info['sequence']

            # 计算相似度
            similarity = self._calculate_similarity(query_seq, target_seq)

            if similarity > best_score:
                best_score = similarity
                best_match = {
                    'cds': seq_info,
                    'locus_tag': locus_tag,
                    'similarity': similarity,
                    'contig': seq_info.get('contig', 'unknown')
                }

        # 只有相似度超过阈值才返回
        if best_match and best_match['similarity'] >= 0.7:
            return best_match

        return None

    def _calculate_similarity(self, seq1, seq2):
        """计算两条蛋白质序列的相似度（优化版：使用 C 加速比对）"""
        if not seq1 or not seq2:
            return 0

        # 使用 Biopython 底层 C 实现的 PairwiseAligner，速度比纯 Python 快得多
        try:
            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'  # 全局比对
            # 设置简单的打分矩阵
            aligner.match_score = 1.0
            aligner.mismatch_score = 0.0
            aligner.open_gap_score = -0.5
            aligner.extend_gap_score = -0.1
            
            score = aligner.score(seq1, seq2)
            max_len = max(len(seq1), len(seq2))
            
            return score / max_len if max_len > 0 else 0
        except Exception:
            # 如果 Biopython 版本过低不支持，回退到简单的长度比
            # 这里做一个极其简单的兜底
            min_len = min(len(seq1), len(seq2))
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
            return matches / min_len if min_len > 0 else 0

    def _match_antismash_batch(self, batch_df, batch_positions):
        """匹配批次结果与antiSMASH（终极版：O(1)极速匹配 + 逻辑修复）"""
        if not self.antismash_cache:
            self.logger.info("   ⚠️ antiSMASH cache不可用，设置N/A")
            # 初始化空列
            cols = ['in_bgc', 'bgc_region_id', 'bgc_type', 'bgc_region_start', 'bgc_region_end',
                    'antismash_gene_start', 'antismash_gene_end', 'antismash_gene_strand',
                    'antismash_locus_tag', 'antismash_product', 'antismash_gene_kind',
                    'matched_contig', 'match_method', 'match_confidence']
            for col in cols:
                batch_df[col] = 'N/A' if col in ['in_bgc', 'bgc_region_id', 'match_method'] else None
            return batch_df

        if batch_df.empty:
            return batch_df

        self.logger.info("   🧬 匹配antiSMASH...")
        self.logger.info("   📍 四步匹配法：ID搜索→坐标验证→坐标搜索→序列比对 (含模糊匹配)")

        # 确定使用的坐标列
        if 'genome_nt_start' in batch_df.columns and batch_df['genome_nt_start'].notna().any():
            start_col = 'genome_nt_start'
            end_col = 'genome_nt_end'
        else:
            start_col = 'start'
            end_col = 'end'

        # 初始化新列
        new_cols = ['in_bgc', 'bgc_region_id', 'bgc_type', 'bgc_region_start', 'bgc_region_end',
                    'antismash_gene_start', 'antismash_gene_end', 'antismash_gene_strand',
                    'antismash_locus_tag', 'antismash_product', 'antismash_gene_kind',
                    'matched_contig', 'match_method', 'match_confidence']
        for col in new_cols:
            batch_df[col] = None
        batch_df['in_bgc'] = 'No'

        # 统计变量
        matched_bgc_count = 0
        match_stats = {'direct_id_match': 0, 'coordinate_match': 0, 'sequence_search': 0, 'total': 0}

        # 辅助函数：标准化名称（移除特殊字符和NZ_前缀）
        def normalize_name(name):
            name = re.sub(r'^NZ_', '', str(name))  # Bug修复：移除NZ_前缀
            return re.sub(r'[^a-zA-Z0-9]', '', name).lower()

        for row in batch_df.itertuples():
            idx = row.Index
            genome = row.genome
            contig = row.contig

            # ✅ 步骤1：获取位置信息
            possible_ids = [f"{genome}___{contig}", contig]
            pos_info = None
            for test_id in possible_ids:
                if test_id in batch_positions:
                    pos_info = batch_positions[test_id]
                    break

            if not pos_info or not pos_info.get('nt_start'):
                batch_df.at[idx, 'match_method'] = 'no_position'
                continue

            # ✅ 步骤2：提取contig名称
            contig_name = pos_info.get('contig')
            if not contig_name:
                contig_name = contig

            try:
                pos_start = int(getattr(row, start_col))
                pos_end = int(getattr(row, end_col))
            except (ValueError, TypeError):
                continue

            if genome not in self.antismash_cache.data:
                continue

            # ✅ 步骤3：Contig 名称匹配 (优化版：使用 Hash Map)
            force_sequence_match = False
            actual_contig_name = None 

            # A. 直接匹配
            if contig_name in self.antismash_cache.data[genome]:
                actual_contig_name = contig_name
            # B. Map 极速查找 (优化点!)
            else:
                norm_target = normalize_name(contig_name)
                # 直接查表，无需循环
                if norm_target in self.antismash_cache.contig_map[genome]:
                    actual_contig_name = self.antismash_cache.contig_map[genome][norm_target]
            
            if not actual_contig_name:
                # 仍然找不到，标记需要强制序列比对
                self.logger.warning(f"   ⚠️ Contig匹配失败: {contig_name}，将尝试全库序列比对")
                force_sequence_match = True
                batch_df.at[idx, 'match_method'] = 'contig_not_found_fallback'
            else:
                # 找到了正确的名字，更新一下
                contig_name = actual_contig_name

            # ✅ 步骤4：执行匹配 (ID -> 坐标 -> 序列)
            matched = False
            
            # 如果 Contig 找到了，尝试 ID 和 坐标匹配
            if not force_sequence_match:
                # A. ID 搜索
                faa_description = pos_info.get('description', '')
                extracted_id = self._extract_id_from_fasta(faa_description)
                
                if extracted_id:
                    converted_locus_tag = self._convert_id_to_locus_tag(extracted_id, contig_name)
                    # 遍历该 Contig 下的所有 CDS
                    cds_list = self.antismash_cache.data[genome][contig_name]['cds']
                    
                    for cds in cds_list:
                        # 尝试匹配 locus_tag
                        if converted_locus_tag and cds['locus_tag'] == converted_locus_tag:
                            # 验证坐标重叠
                            overlap_start = max(pos_start, cds['start'])
                            overlap_end = min(pos_end, cds['end'])
                            if overlap_end - overlap_start > 0:
                                matched = True
                                batch_df.at[idx, 'match_method'] = 'direct_id_match'
                                batch_df.at[idx, 'match_confidence'] = 'high'
                                self._record_match(batch_df, idx, cds, contig_name)
                                match_stats['direct_id_match'] += 1
                                break
                
                # B. 坐标搜索 (如果 ID 没匹配上)
                if not matched:
                    overlapping_cds = self.antismash_cache.find_overlapping_cds(
                        genome, contig_name, pos_start, pos_end
                    )
                    # 找重叠最大的
                    best_cds = None
                    max_overlap = 0
                    for cds in overlapping_cds:
                        ov_len = min(pos_end, cds['end']) - max(pos_start, cds['start'])
                        if ov_len > max_overlap:
                            max_overlap = ov_len
                            best_cds = cds
                    
                    if best_cds and max_overlap > 0:
                        matched = True
                        batch_df.at[idx, 'match_method'] = 'coordinate_match'
                        batch_df.at[idx, 'match_confidence'] = 'medium'
                        self._record_match(batch_df, idx, best_cds, contig_name)
                        match_stats['coordinate_match'] += 1

            # C. 序列比对 (如果前面的都没匹配上，或者 Contig 根本没找到)
            if not matched or force_sequence_match:
                query_seq = self._get_query_sequence(row.query_id)
                if query_seq:
                    # 如果是 force_sequence_match，则在全库搜索；否则只在当前 Contig 搜
                    gbk_sequences = self._load_gbk_sequences(
                        genome, contig_name if not force_sequence_match else None, 
                        force_all_contigs=force_sequence_match
                    )
                    
                    best_match = self._perform_sequence_alignment(query_seq, gbk_sequences)
                    if best_match:
                        matched = True
                        batch_df.at[idx, 'match_method'] = 'sequence_search'
                        batch_df.at[idx, 'match_confidence'] = 'medium' if force_sequence_match else 'low'
                        # 注意：这里会更新 contig_name 为序列比对找到的那个 contig
                        self._record_match(batch_df, idx, best_match['cds'], best_match['contig'])
                        match_stats['sequence_search'] += 1
                        
                        # 【重要修复】强制更新 contig_name 变量，供步骤5使用
                        if force_sequence_match:
                            contig_name = best_match['contig']

            # ✅ 步骤5：检查是否在 BGC Region 中 (逻辑修复)
            # 优先使用 matched_contig (如果有)，否则使用当前的 contig_name
            final_contig = batch_df.at[idx, 'matched_contig']
            if not final_contig:
                final_contig = contig_name

            # 只有当 contig 名称存在于缓存中时才检查
            if final_contig and final_contig in self.antismash_cache.data[genome]:
                data = self.antismash_cache.data[genome][final_contig]
                tolerance = 100
                in_region = False
                
                # 如果没有精确匹配到CDS，使用原始搜索位置；如果匹配到了CDS，使用CDS位置
                chk_start = batch_df.at[idx, 'antismash_gene_start'] if matched else pos_start
                chk_end = batch_df.at[idx, 'antismash_gene_end'] if matched else pos_end
                
                if chk_start and chk_end:
                    for region in data['regions']:
                        if (chk_start - tolerance <= region['end'] and 
                            chk_end + tolerance >= region['start']):
                            in_region = True
                            matched_bgc_count += 1
                            batch_df.at[idx, 'in_bgc'] = 'Yes'
                            batch_df.at[idx, 'bgc_region_id'] = f"region{region['region_number']}" if region['region_number'] else f"region_{region['start']}_{region['end']}"
                            batch_df.at[idx, 'bgc_type'] = region['type']
                            batch_df.at[idx, 'bgc_region_start'] = region['start']
                            batch_df.at[idx, 'bgc_region_end'] = region['end']
                            if pd.isna(batch_df.at[idx, 'matched_contig']):
                                batch_df.at[idx, 'matched_contig'] = final_contig
                            break
                
                if not in_region:
                     batch_df.at[idx, 'in_bgc'] = 'No'
            else:
                 batch_df.at[idx, 'in_bgc'] = 'No (contig_not_found)'

            match_stats['total'] += 1

        self.logger.info(f"   ✓ 在BGC内: {matched_bgc_count} 条")
        return batch_df

    def _record_match(self, df, idx, cds, contig):
        """辅助函数：记录匹配信息到DataFrame"""
        df.at[idx, 'antismash_gene_start'] = cds['start']
        df.at[idx, 'antismash_gene_end'] = cds['end']
        df.at[idx, 'antismash_gene_strand'] = cds['strand']
        df.at[idx, 'antismash_locus_tag'] = cds['locus_tag']
        df.at[idx, 'antismash_product'] = cds['product']
        df.at[idx, 'antismash_gene_kind'] = cds.get('gene_kind', '')
        df.at[idx, 'matched_contig'] = contig

    def run(self):
        """运行完整流程 - 流式处理版本"""
        self.logger.info("\n" + "="*80)
        self.logger.info("🚀 开始流式处理流程（每批次完整处理）")
        self.logger.info("="*80)
        
        total_start_time = time.time()
        batch_output_files = []  # 记录每个批次的输出文件
        
        try:
            for batch in self.batch_plan:
                # === 新增：断点续传 ===
                batch_output_file = self.output_dir / "batches" / f"batch_{batch['batch_id']}_complete.tsv"
                # 如果文件存在且不为空，说明这个批次以前跑成功过，直接跳过
                if batch_output_file.exists() and batch_output_file.stat().st_size > 100:
                    self.logger.info(f"⏭️ 批次 {batch['batch_id']} 已完成，跳过...")
                    batch_output_files.append(batch_output_file)
                    continue
                # ====================

                batch_start = time.time()
                self.logger.info(f"\n{'='*20} 批次 {batch['batch_id']}/{len(self.batch_plan)} {'='*20}")
                self.logger.info(f"💾 当前内存: {get_memory_usage():.1f} MB")

                try:
                    # 1️⃣ 搜索
                    self.logger.info("1️⃣ 搜索批次...")
                    db_path, genome_mapping, batch_seq_file = self._create_batch_database(
                        batch['items'], batch['batch_id']
                    )

                    batch_results = self._search_batch(db_path, genome_mapping, batch['batch_id'])

                    if not batch_results:
                        self.logger.info("   ⚠️ 该批次无匹配结果")
                        self._cleanup_batch_files(batch_seq_file, db_path, batch['batch_id'])
                        continue

                    batch_df = pd.DataFrame(batch_results)
                    self.logger.info(f"   找到 {len(batch_df)} 个匹配")

                    # 过滤
                    original_count = len(batch_df)
                    batch_df = batch_df[
                        (batch_df['identity'] >= self.min_identity) &
                        (batch_df['coverage'] >= self.min_coverage)
                    ]
                    self.logger.info(f"   过滤后: {len(batch_df)} 条 (Identity≥{self.min_identity}%, Coverage≥{self.min_coverage}%)")

                    if batch_df.empty:
                        self.logger.info("   ⚠️ 过滤后无结果")
                        self._cleanup_batch_files(batch_seq_file, db_path, batch['batch_id'])
                        continue

                    # 2️⃣ 位置修复（加载批次位置信息）
                    self.logger.info("2️⃣ 处理批次...")
                    try:
                        batch_positions = self._load_batch_positions(batch['items'])
                        batch_df = self._fix_positions_batch(batch_df, batch_positions)
                    except Exception as e:
                        self.logger.error(f"   ⚠️ 位置修复失败: {e}")
                        self.logger.info(f"   ℹ️ 跳过位置修复，继续处理批次 {batch['batch_id']}")
                        # 设置默认位置信息
                        batch_df['genome_nt_start'] = None
                        batch_df['genome_nt_end'] = None
                        batch_df['gene_nt_location'] = None
                        batch_df['gene_strand'] = None
                        batch_df['position_fixed'] = False

                    # 3️⃣ antiSMASH匹配（加载批次antiSMASH数据）
                    try:
                        self.antismash_cache = self._load_batch_antismash(batch['items'])
                        batch_df = self._match_antismash_batch(batch_df, batch_positions)
                    except Exception as e:
                        self.logger.error(f"   ⚠️ antiSMASH处理失败: {e}")
                        self.logger.info(f"   ℹ️ 跳过antiSMASH匹配，继续处理批次 {batch['batch_id']}")
                        # 设置默认antiSMASH列
                        batch_df['in_bgc'] = 'N/A'
                        batch_df['bgc_region_id'] = 'N/A'
                        batch_df['bgc_type'] = 'N/A'
                        batch_df['bgc_region_start'] = None
                        batch_df['bgc_region_end'] = None
                        batch_df['antismash_gene_start'] = None
                        batch_df['antismash_gene_end'] = None
                        batch_df['antismash_gene_strand'] = None
                        batch_df['antismash_locus_tag'] = None
                        batch_df['antismash_product'] = None
                        batch_df['match_method'] = 'N/A'
                        batch_df['match_confidence'] = 'N/A'

                    # 4️⃣ 保存批次结果
                    try:
                        # 确保 batch_output_file 变量已定义（使用之前用于判断是否存在的路径）
                        batch_df.to_csv(batch_output_file, sep='\t', index=False)
                        batch_output_files.append(batch_output_file)
                        self.logger.info(f"   ✓ 保存: {batch_output_file.name}")
                    except Exception as e:
                        self.logger.error(f"   ⚠️ 保存批次结果失败: {e}")
                        self.logger.info(f"   ℹ️ 跳过批次 {batch['batch_id']} 的结果保存")

                except Exception as e:
                    self.logger.error(f"   ❌ 批次 {batch['batch_id']} 处理失败: {e}")
                    self.logger.info(f"   ℹ️ 跳过批次 {batch['batch_id']}")
                    # 继续处理下一个批次
                    continue

                finally:
                    # 5️⃣ 清理临时文件和内存（无论成功失败都执行）
                    try:
                        if 'batch_seq_file' in locals():
                            self._cleanup_batch_files(batch_seq_file, db_path, batch['batch_id'])
                        if 'batch_results' in locals():
                            del batch_results
                        if 'batch_df' in locals():
                            del batch_df
                        if 'genome_mapping' in locals():
                            del genome_mapping
                        if 'batch_positions' in locals():
                            del batch_positions
                        self.antismash_cache = None  # 释放批次的antiSMASH缓存
                        gc.collect()
                    except Exception as e:
                        self.logger.warning(f"   ⚠️ 清理临时数据时出错: {e}")

                batch_time = time.time() - batch_start
                self.logger.info(f"   ⏱️ 批次耗时: {batch_time:.1f} 秒")
                self.logger.info(f"   💾 当前内存: {get_memory_usage():.1f} MB")
            
            # 6️⃣ 合并所有批次结果
            self.logger.info("\n" + "="*80)
            self.logger.info("📊 合并所有批次结果")
            self.logger.info("="*80)
            
            if not batch_output_files:
                self.logger.warning("没有任何批次产生结果")
                return False
            
            final_results = self._merge_batch_results(batch_output_files)
            
            # 7️⃣ 生成报告
            total_time = time.time() - total_start_time
            self._generate_detailed_report(final_results, total_time)
            
            # 8️⃣ 打印总结
            self._print_summary(final_results, total_time)
            
            return True
            
        except MemoryError as e:
            self.logger.error(f"内存错误: {e}")
            self.logger.error("建议措施:")
            self.logger.error("  1. 减少 --batch-size 参数")
            self.logger.error("  2. 增加系统内存或swap空间")
            return False
        except Exception as e:
            self.logger.error(f"流程失败: {e}")
            import traceback
            traceback.print_exc()
            return False
        
        finally:
            # 清理临时目录
            if self.temp_dir and self.temp_dir.exists():
                try:
                    shutil.rmtree(self.temp_dir)
                    self.logger.info(f"✓ 清理临时目录: {self.temp_dir}")
                except Exception as e:
                    self.logger.warning(f"清理临时目录失败: {e}")
    
    def _merge_batch_results(self, batch_files):
        """合并批次结果（文件级操作，内存友好）"""
        self.logger.info(f"合并 {len(batch_files)} 个批次文件...")
        
        # 方法1：如果文件不太大，直接concat
        total_size = sum(f.stat().st_size for f in batch_files) / (1024 * 1024)
        self.logger.info(f"总文件大小: {total_size:.1f} MB")
        
        if total_size < 500:  # 小于500MB，直接读取
            dfs = []
            for batch_file in batch_files:
                df = pd.read_csv(batch_file, sep='\t')
                dfs.append(df)
                self.logger.info(f"  读取 {batch_file.name}: {len(df)} 条")
            
            final_df = pd.concat(dfs, ignore_index=True)
            del dfs
            gc.collect()
        
        else:  # 大文件，使用文件拼接
            self.logger.info("使用文件拼接模式...")
            final_file = self.output_dir / "reports" / "complete_results.tsv"
            
            with open(final_file, 'w') as outf:
                for i, batch_file in enumerate(batch_files):
                    with open(batch_file, 'r') as inf:
                        if i == 0:
                            # 第一个文件：包含header
                            outf.write(inf.read())
                        else:
                            # 后续文件：跳过header
                            next(inf)  # 跳过header行
                            outf.write(inf.read())
            
            # 读取合并后的结果
            final_df = pd.read_csv(final_file, sep='\t')
        
        self.logger.info(f"✓ 合并完成: {len(final_df)} 条记录")

        # 去重：基于位置（genome + contig + start + end）
        # 因为参考序列中有相似序列可能匹配到同一位置
        before_dedup = len(final_df)

        # 确定位置列名
        pos_cols = []
        if 'genome' in final_df.columns:
            pos_cols.append('genome')
        if 'matched_contig' in final_df.columns:
            pos_cols.append('matched_contig')
        elif 'contig' in final_df.columns:
            pos_cols.append('contig')

        # 使用核苷酸位置去重（更准确）
        if 'genome_nt_start' in final_df.columns and 'genome_nt_end' in final_df.columns:
            pos_cols.extend(['genome_nt_start', 'genome_nt_end'])
        elif 'start' in final_df.columns and 'end' in final_df.columns:
            pos_cols.extend(['start', 'end'])

        if len(pos_cols) >= 3:  # 至少有genome + contig + start
            # 保留每个位置的最佳匹配（按identity或evalue排序）
            if 'identity' in final_df.columns:
                final_df = final_df.sort_values('identity', ascending=False)
            elif 'evalue' in final_df.columns:
                final_df = final_df.sort_values('evalue', ascending=True)

            final_df = final_df.drop_duplicates(subset=pos_cols, keep='first')
            after_dedup = len(final_df)

            if before_dedup > after_dedup:
                self.logger.info(f"✓ 位置去重: {before_dedup} → {after_dedup} 条 (移除 {before_dedup - after_dedup} 重复)")

        return final_df
    
    def _generate_detailed_report(self, results_df, total_time):
        """生成详细报告"""
        self.logger.info("\n📊 生成详细报告...")

        # 保存完整结果
        results_file = self.output_dir / "reports" / "complete_results.tsv"
        results_df.to_csv(results_file, sep='\t', index=False)
        self.logger.info(f"✓ 完整结果: {results_file}")

        # 保存BGC内的结果
        if 'in_bgc' in results_df.columns:
            bgc_hits = results_df[results_df['in_bgc'] == 'Yes']
            if not bgc_hits.empty:
                bgc_file = self.output_dir / "reports" / "bgc_hits.tsv"
                bgc_hits.to_csv(bgc_file, sep='\t', index=False)
                self.logger.info(f"✓ BGC内结果: {bgc_file} ({len(bgc_hits)} 条)")

        # 保存按匹配方法分类的结果
        if 'match_method' in results_df.columns:
            for method in results_df['match_method'].unique():
                if pd.notna(method) and method not in ['N/A', None]:
                    method_hits = results_df[results_df['match_method'] == method]
                    method_file = self.output_dir / "reports" / f"match_method_{method}.tsv"
                    method_hits.to_csv(method_file, sep='\t', index=False)
                    self.logger.info(f"✓ {method}结果: {method_file} ({len(method_hits)} 条)")

        # 生成Markdown报告（增强版）
        report_content = f"""# 流式处理搜索报告

## 📊 搜索概要
- **搜索时间**: {time.strftime('%Y-%m-%d %H:%M:%S')}
- **总耗时**: {total_time:.1f} 秒
- **查询类型**: {self.query_type.upper()}
- **数据库类型**: {self.database_type.upper()}
- **搜索策略**: {self.blast_mode.upper()}
- **总结果数**: {len(results_df)} 条

## 🔧 处理模式
- ✅ **流式处理**: 每批次完整处理后释放内存
- ✅ **批次数量**: {len(self.batch_plan)}
- ✅ **批次大小**: ~{self.batch_size} 个文件/批次
- ✅ **增强ID匹配**: 四步匹配法（ID搜索→坐标验证→坐标搜索→序列比对）

## 📈 统计
"""

        if self.is_protein_search and 'position_fixed' in results_df.columns:
            fixed_count = results_df['position_fixed'].sum()
            report_content += f"""
### 位置修复
- 成功转换: {fixed_count}/{len(results_df)} 条 ({fixed_count/len(results_df)*100:.1f}%)
"""

        if 'in_bgc' in results_df.columns:
            in_bgc_count = len(results_df[results_df['in_bgc'] == 'Yes'])
            report_content += f"""
### antiSMASH匹配
- 位于BGC内: {in_bgc_count}/{len(results_df)} 条 ({in_bgc_count/len(results_df)*100:.1f}%)
"""

        # 添加匹配方法统计
        if 'match_method' in results_df.columns:
            match_counts = results_df['match_method'].value_counts()
            report_content += f"""
### 匹配方法统计
"""

            for method, count in match_counts.items():
                if pd.notna(method):
                    percentage = count / len(results_df) * 100
                    description = ''
                    if method == 'direct_id_match':
                        description = 'ID搜索+坐标验证'
                    elif method == 'coordinate_match':
                        description = '坐标搜索'
                    elif method == 'sequence_search':
                        description = '序列比对匹配'
                    elif method == 'contig_not_found_sequence_fallback':
                        description = 'Contig匹配失败后的序列比对兜底'

                    report_content += f"- **{method}**: {count} 条 ({percentage:.1f}%) - {description}\n"

            # 添加匹配置信度统计
            if 'match_confidence' in results_df.columns:
                confidence_counts = results_df['match_confidence'].value_counts()
                report_content += f"""
### 匹配置信度分布
"""
                for confidence, count in confidence_counts.items():
                    if pd.notna(confidence):
                        percentage = count / len(results_df) * 100
                        report_content += f"- **{confidence}**: {count} 条 ({percentage:.1f}%)\n"

        report_file = self.output_dir / "reports" / "detailed_report.md"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report_content)
        self.logger.info(f"✓ 详细报告: {report_file}")
    
    def _print_summary(self, results_df, total_time):
        """打印最终总结"""
        print("\n" + "="*80)
        print("🎉 流式处理完成！")
        print("="*80)
        print(f"⏱️ 总耗时: {total_time:.1f} 秒")
        print(f"📊 处理模式: 流式处理 ({len(self.batch_plan)} 个批次)")
        print(f"✅ 总结果数: {len(results_df)} 条")
        
        if self.is_protein_search and 'position_fixed' in results_df.columns:
            fixed = results_df['position_fixed'].sum()
            print(f"🔧 位置修复: {fixed}/{len(results_df)} 条")
        
        if 'in_bgc' in results_df.columns:
            in_bgc = len(results_df[results_df['in_bgc'] == 'Yes'])
            print(f"🧬 BGC内命中: {in_bgc}/{len(results_df)} 条")
        
        print(f"\n📁 输出文件:")
        print(f"   - 完整结果: reports/complete_results.tsv")
        print(f"   - 详细报告: reports/detailed_report.md")
        
        if 'in_bgc' in results_df.columns and len(results_df[results_df['in_bgc'] == 'Yes']) > 0:
            print(f"   - BGC命中: reports/bgc_hits.tsv")
        
        print("="*80)

# ============================================================================
# 主函数
# ============================================================================

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="流式处理版本：整合的序列搜索和BGC匹配完整流程（使用整个基因位置匹配）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
核心改进：
  ✅ 每批次走完整流程（搜索→修复→匹配→保存→释放内存）
  ✅ sequence_positions和antismash_cache预加载一次，全程复用
  ✅ 批次结果立即处理并保存，不在内存中累积
  ✅ 大幅降低内存峰值，避免内存爆炸
  ✅ 使用整个基因位置匹配BGC（更简单、更准确）

使用示例:
  # 基本搜索（不含antiSMASH）
  python conmon_search_pipeline_stream.py \\
    -d /path/to/database/ \\
    -q /path/to/queries/ \\
    --batch-size 50
  
  # 完整流程（含antiSMASH匹配）
  python conmon_search_pipeline_stream.py \\
    -d /path/to/database/ \\
    -q /path/to/queries/ \\
    -a /path/to/antismash_results/ \\
    --batch-size 50 \\
    --threads 8
        """
    )
    
    parser.add_argument("-d", "--database", required=True,
                       help="数据库序列文件目录")
    parser.add_argument("-q", "--queries", required=True,
                       help="查询序列文件目录")
    parser.add_argument("-a", "--antismash",
                       help="antiSMASH结果目录（可选）")
    parser.add_argument("-o", "--output", default="search_results",
                       help="输出目录（默认: search_results）")
    parser.add_argument("-b", "--batch-size", type=int, default=100,
                       help="每批处理的文件数量（默认: 50，建议范围: 20-100）")
    parser.add_argument("-t", "--threads", type=int, default=4,
                       help="并行线程数（默认: 4）")
    parser.add_argument("--use-diamond", action='store_true', default=True,
                       help="使用DIAMOND搜索（默认）")
    parser.add_argument("--use-blast", action='store_true',
                       help="使用BLAST搜索")
    parser.add_argument("--blast-only", action='store_true',
                       help="仅使用BLAST（禁用DIAMOND）")
    parser.add_argument("--min-identity", type=float, default=50,
                       help="最小相似度阈值（默认: 50%%）")
    parser.add_argument("--min-coverage", type=float, default=50,
                       help="最小覆盖度阈值（默认: 50%%）")
    
    args = parser.parse_args()
    
    # 验证输入
    if not Path(args.database).exists():
        print(f"❌ 数据库目录不存在: {args.database}")
        sys.exit(1)
    
    if not Path(args.queries).exists():
        print(f"❌ 查询目录不存在: {args.queries}")
        sys.exit(1)
    
    # antiSMASH目录验证
    antismash_path = None
    if args.antismash:
        antismash_path = Path(args.antismash)
        if not antismash_path.exists():
            print(f"⚠️ antiSMASH目录不存在: {args.antismash}")
            print("   将跳过BGC匹配")
            antismash_path = None
        else:
            print(f"✓ antiSMASH目录: {antismash_path}")
    
    # 确定搜索方法
    use_diamond = args.use_diamond and not args.blast_only
    use_blast = args.use_blast or args.blast_only
    
    if not use_diamond and not use_blast:
        print("⚠️ 未指定搜索方法，默认使用DIAMOND")
        use_diamond = True
    
    # 创建流程对象并运行
    pipeline = StreamProcessingPipeline(
        database_dir=args.database,
        query_dir=args.queries,
        output_dir=args.output,
        batch_size=args.batch_size,
        threads=args.threads,
        use_diamond=use_diamond,
        use_blast=use_blast,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        antismash_dir=antismash_path
    )
    
    success = pipeline.run()
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()