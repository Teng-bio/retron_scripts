#!/usr/bin/env python3
"""
Retron系统挖掘流程 v2.0 - 主控脚本

基于Mestre et al., 2020 NAR的方法论，采用统计关联分析策略：
1. RT搜索 (DIAMOND/BLAST)
2. Motif过滤 (NAXXH + VTG)
3. RT系统发育分析 (排除Group II/DGR/CRISPR-RT)
4. 邻近蛋白提取 (±30kb, Mestre方法)
5. MMseqs2蛋白聚类
6. Phyvalue统计关联分析
7. 综合报告与Retron类型分类
8. ncRNA共变验证 (可选，多序列比对方法)

方法参考：
- Mestre et al., 2020 NAR: 统计关联分析，基于Phyvalue筛选显著关联蛋白簇
- 与Millman方法的区别：本流程不依赖DefenseFinder，而是通过无监督聚类发现关联蛋白
- ncRNA验证采用共变分析，适用于高GC生物
"""

import argparse
import subprocess
import sys
import yaml
from pathlib import Path
import logging
import shutil
from datetime import datetime


def setup_logging(output_dir, verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    log_file = output_dir / "logs" / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def load_config(config_file):
    """加载配置文件"""
    with open(config_file, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def run_step(cmd, step_name, logger, dry_run=False):
    """运行单个步骤"""
    logger.info(f"\n{'='*60}")
    logger.info(f"▶ {step_name}")
    logger.info(f"{'='*60}")
    logger.debug(f"命令: {' '.join(cmd)}")

    if dry_run:
        logger.info("[DRY RUN] 跳过执行")
        return True

    try:
        result = subprocess.run(cmd, check=True)
        logger.info(f"✓ {step_name} 完成")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"✗ {step_name} 失败: {e}")
        return False
    except FileNotFoundError:
        logger.error(f"✗ 脚本不存在: {cmd[0]}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Retron系统挖掘流程 v2.0 (Mestre方法)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 使用配置文件
  python main_pipeline.py -c config.yaml

  # 预览命令
  python main_pipeline.py -c config.yaml --dry-run

  # 命令行参数
  python main_pipeline.py \\
    --search-results /path/to/complete_results.tsv \\
    --faa-dir /path/to/faa \\
    --fasta-dir /path/to/fasta \\
    --antismash-dir /path/to/antismash \\
    --output ./retron_v2_results

流程概述:
  Step 1: RT搜索 (DIAMOND)
  Step 2: Motif过滤 (NAXXH+VTG)
  Step 3: RT系统发育分析 (MAFFT+FastTree, 排除非Retron RT)
  Step 4: 邻近蛋白提取 (±30kb)
  Step 5: MMseqs2蛋白聚类
  Step 6: Phyvalue统计关联分析
  Step 7: 综合报告与类型分类
  Step 8: ncRNA共变验证 (可选, MAFFT+RNAalifold+R-scape)
        """
    )

    # 配置文件
    parser.add_argument('-c', '--config', help='YAML配置文件')

    # 输入路径
    parser.add_argument('--query-dir', help='RT查询序列目录 (Step 1)')
    parser.add_argument('--search-results', help='已有的搜索结果 (跳过Step 1)')
    parser.add_argument('--faa-dir', help='蛋白序列目录 (.faa)')
    parser.add_argument('--fasta-dir', help='基因组序列目录 (.fasta)')
    parser.add_argument('--antismash-dir', help='antiSMASH结果目录')

    # 输出
    parser.add_argument('-o', '--output',
                        help='输出目录 (默认从配置文件读取)')

    # 流程控制
    parser.add_argument('--start-step', type=int, default=1, help='起始步骤 (1-8)')
    parser.add_argument('--end-step', type=int, default=8, help='结束步骤 (1-8)')
    parser.add_argument('--skip-ncrna', action='store_true',
                        help='跳过Step 8 ncRNA共变验证')

    # 搜索参数 (Step 1)
    parser.add_argument('--min-identity', type=float, default=25,
                        help='最小序列相似度%% (默认: 25)')
    parser.add_argument('--min-coverage', type=float, default=40,
                        help='最小覆盖度%% (默认: 40)')

    # 系统发育参数 (Step 3)
    parser.add_argument('--distance-threshold', type=float, default=0.5,
                        help='系统发育距离阈值 (默认: 0.5)')

    # 邻近蛋白提取参数 (Step 4)
    parser.add_argument('--neighbor-distance', type=int, default=30,
                        help='邻近蛋白提取距离，单位kb (默认: 30)')

    # MMseqs2 聚类参数 (Step 5)
    parser.add_argument('--min-seq-id', type=float, default=0.3,
                        help='MMseqs2最小序列一致性 (默认: 0.3 = 30%%)')
    parser.add_argument('--coverage', type=float, default=0.8,
                        help='MMseqs2最小覆盖度 (默认: 0.8 = 80%%)')

    # Phyvalue 统计参数 (Step 6)
    parser.add_argument('--n-permutations', type=int, default=1000,
                        help='置换检验次数 (默认: 1000)')
    parser.add_argument('--min-occurrence', type=int, default=5,
                        help='最小出现次数阈值 (默认: 5)')
    parser.add_argument('--min-phyvalue', type=float, default=2.0,
                        help='最小Phyvalue阈值 (默认: 2.0)')

    # ncRNA验证参数 (Step 8)
    parser.add_argument('--upstream', type=int, default=600,
                        help='ncRNA上游提取长度 (默认: 600bp)')
    parser.add_argument('--min-group-size', type=int, default=3,
                        help='ncRNA分析最小分组大小 (默认: 3)')

    # 通用参数
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--dry-run', action='store_true', help='只显示命令')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    # 加载配置
    config = {}
    if args.config:
        config_path = Path(args.config)
        if config_path.exists():
            config = load_config(config_path)
            print(f"✓ 已加载配置文件: {config_path.absolute()}")
        else:
            print(f"⚠️ 配置文件不存在: {config_path.absolute()}")
            print(f"   将使用默认配置")

    # 合并配置
    def get_config(key, default=None):
        """从配置或参数获取值"""
        keys = key.split('.')
        val = config
        for k in keys:
            if isinstance(val, dict):
                val = val.get(k)
            else:
                val = None
                break
        return val if val is not None else default

    # 路径配置
    query_dir = args.query_dir or get_config('paths.query_dir')
    search_results = args.search_results or get_config('pipeline.existing_search_results')
    faa_dir = args.faa_dir or get_config('paths.faa_dir')
    fasta_dir = args.fasta_dir or get_config('paths.fasta_dir')
    antismash_dir = args.antismash_dir or get_config('paths.antismash_dir')
    reference_dir = get_config('paths.reference_dir', '/home/teng/claude_code/retron/database/rt_reference')
    output_dir = Path(args.output) if args.output else Path(get_config('paths.output_dir', './retron_v2_results'))

    # 参数配置 (命令行优先级高于配置文件)
    min_identity = args.min_identity if args.min_identity != 25 else get_config('search.min_identity', 25)
    min_coverage = args.min_coverage if args.min_coverage != 40 else get_config('search.min_coverage', 40)
    threads = args.threads if args.threads != 8 else get_config('general.threads', 8)

    # 系统发育参数
    distance_threshold = args.distance_threshold if args.distance_threshold != 0.5 else get_config('phylogeny.distance_threshold', 0.5)

    # 邻近蛋白参数
    neighbor_distance = args.neighbor_distance if args.neighbor_distance != 30 else get_config('neighbors.distance', 30)

    # 聚类参数
    min_seq_id = args.min_seq_id if args.min_seq_id != 0.3 else get_config('clustering.min_seq_id', 0.3)
    coverage = args.coverage if args.coverage != 0.8 else get_config('clustering.coverage', 0.8)

    # 关联分析参数
    n_permutations = args.n_permutations if args.n_permutations != 1000 else get_config('association.n_permutations', 1000)
    min_occurrence = args.min_occurrence if args.min_occurrence != 5 else get_config('association.min_occurrence', 5)
    min_phyvalue = args.min_phyvalue if args.min_phyvalue != 2.0 else get_config('association.min_phyvalue', 2.0)

    # ncRNA参数
    upstream = args.upstream if args.upstream != 600 else get_config('ncrna.upstream', 600)
    min_group_size = args.min_group_size if args.min_group_size != 3 else get_config('ncrna.min_group_size', 3)

    start_step = args.start_step or get_config('pipeline.start_step', 1)
    end_step = args.end_step or get_config('pipeline.end_step', 8)

    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging(output_dir, args.verbose)

    # 脚本目录
    script_dir = Path(__file__).parent

    # 步骤输出目录
    dirs = {
        1: output_dir / "01_search",
        2: output_dir / "02_motif_filtered",
        3: output_dir / "03_phylogeny",
        4: output_dir / "04_neighbors",
        5: output_dir / "05_clustering",
        6: output_dir / "06_association",
        7: output_dir / "07_report",
        8: output_dir / "08_ncrna"
    }

    logger.info("=" * 60)
    logger.info("Retron系统挖掘流程 v2.0 (Mestre方法)")
    logger.info("=" * 60)
    logger.info(f"输出目录: {output_dir}")
    logger.info(f"步骤范围: {start_step} - {end_step}")

    # 显示配置信息
    if config:
        logger.info(f"配置文件: {args.config}")

    success = True

    # 步骤输出文件映射
    step_outputs = {
        1: dirs[1] / "reports" / "complete_results.tsv",
        2: dirs[2] / "motif_filtered.tsv",
        3: dirs[3] / "retron_candidates.tsv",
        4: dirs[4] / "rt_info.tsv",
        5: dirs[5] / "rt_cluster_matrix.tsv",
        6: dirs[6] / "phyvalue_analysis.tsv",
        7: dirs[7] / "final_candidates.tsv",
    }

    # 确定初始输入文件
    if start_step > 1:
        prev_step = start_step - 1
        if prev_step in step_outputs and step_outputs[prev_step].exists():
            current_input = step_outputs[prev_step]
            logger.info(f"从Step {start_step}开始，使用输入: {current_input}")
        else:
            current_input = search_results if search_results else step_outputs.get(1)
    else:
        current_input = search_results

    # ==================== Step 1: RT搜索 ====================
    if start_step <= 1 <= end_step and not search_results:
        if not query_dir or not faa_dir:
            logger.error("Step 1 需要 --query-dir 和 --faa-dir")
            sys.exit(1)

        cmd = [
            sys.executable, str(script_dir / "step01_search.py"),
            "-q", str(query_dir),
            "-d", str(faa_dir),
            "-o", str(dirs[1]),
            "-t", str(threads),
            "--min-identity", str(min_identity),
            "--min-coverage", str(min_coverage)
        ]
        if antismash_dir:
            cmd.extend(["-a", str(antismash_dir)])

        success = run_step(cmd, "Step 1: RT搜索", logger, args.dry_run) and success
        current_input = dirs[1] / "reports" / "complete_results.tsv"
    elif start_step <= 1:
        current_input = Path(search_results) if search_results else dirs[1] / "reports" / "complete_results.tsv"

    # ==================== Step 2: Motif过滤 ====================
    if start_step <= 2 <= end_step and success:
        if not faa_dir:
            logger.error("Step 2 需要 --faa-dir")
            sys.exit(1)

        cmd = [
            sys.executable, str(script_dir / "step02_motif_filter.py"),
            "-i", str(current_input),
            "-f", str(faa_dir),
            "-o", str(dirs[2])
        ]
        success = run_step(cmd, "Step 2: Motif过滤 (NAXXH/VTG)", logger, args.dry_run) and success
        current_input = dirs[2] / "motif_filtered.tsv"

    # ==================== Step 3: RT系统发育分析 ====================
    if start_step <= 3 <= end_step and success:
        if not faa_dir:
            logger.error("Step 3 需要 --faa-dir")
            sys.exit(1)

        # 检查工具
        if not shutil.which('mafft'):
            logger.warning("MAFFT未安装，请先安装: conda install -c bioconda mafft")
        if not shutil.which('fasttree'):
            logger.warning("FastTree未安装，请先安装: conda install -c bioconda fasttree")

        cmd = [
            sys.executable, str(script_dir / "step03_phylogeny.py"),
            "-i", str(current_input),
            "-f", str(faa_dir),
            "-r", str(reference_dir),
            "-o", str(dirs[3]),
            "--distance-threshold", str(distance_threshold),
            "--threads", str(threads)
        ]
        if args.verbose:
            cmd.append("-v")

        success = run_step(cmd, "Step 3: RT系统发育分析", logger, args.dry_run) and success
        current_input = dirs[3] / "retron_candidates.tsv"

    # ==================== Step 4: 邻近蛋白提取 ====================
    if start_step <= 4 <= end_step and success:
        if not antismash_dir:
            logger.error("Step 4 需要 --antismash-dir")
            sys.exit(1)

        cmd = [
            sys.executable, str(script_dir / "step04_extract_neighbors.py"),
            "-i", str(current_input),
            "-a", str(antismash_dir),
            "-o", str(dirs[4]),
            "--distance", str(neighbor_distance)
        ]
        if args.verbose:
            cmd.append("-v")

        success = run_step(cmd, f"Step 4: 邻近蛋白提取 (±{neighbor_distance}kb)", logger, args.dry_run) and success

    # ==================== Step 5: MMseqs2蛋白聚类 ====================
    if start_step <= 5 <= end_step and success:
        step04_dir = dirs[4]
        if not (step04_dir / "all_neighbors.faa").exists() and not args.dry_run:
            logger.warning("Step 4输出不存在，跳过Step 5")
        else:
            if not shutil.which('mmseqs'):
                logger.warning("MMseqs2未安装，请先安装: conda install -c bioconda mmseqs2")

            cmd = [
                sys.executable, str(script_dir / "step05_mmseqs_cluster.py"),
                "-i", str(step04_dir),
                "-o", str(dirs[5]),
                "--min-seq-id", str(min_seq_id),
                "--coverage", str(coverage),
                "--threads", str(threads)
            ]
            if args.verbose:
                cmd.append("-v")

            success = run_step(cmd, "Step 5: MMseqs2蛋白聚类", logger, args.dry_run) and success

    # ==================== Step 6: Phyvalue统计关联分析 ====================
    if start_step <= 6 <= end_step and success:
        step05_dir = dirs[5]
        if not (step05_dir / "rt_cluster_matrix.tsv").exists() and not args.dry_run:
            logger.warning("Step 5输出不存在，跳过Step 6")
        else:
            cmd = [
                sys.executable, str(script_dir / "step06_phyvalue_analysis.py"),
                "-i", str(step05_dir),
                "-o", str(dirs[6]),
                "--n-permutations", str(n_permutations),
                "--min-occurrence", str(min_occurrence),
                "--min-phyvalue", str(min_phyvalue)
            ]
            if args.verbose:
                cmd.append("-v")

            success = run_step(cmd, "Step 6: Phyvalue统计关联分析", logger, args.dry_run) and success

    # ==================== Step 7: 综合报告与类型分类 ====================
    if start_step <= 7 <= end_step and success:
        cmd = [
            sys.executable, str(script_dir / "step07_final_report.py"),
            "-o", str(dirs[7]),
        ]

        # 添加可选输入文件
        optional_inputs = [
            ("--rt-info", dirs[4] / "rt_info.tsv"),
            ("--motif", dirs[2] / "motif_filtered.tsv"),
            ("--phylogeny", dirs[3] / "retron_candidates.tsv"),
            ("--cluster-stats", dirs[6] / "phyvalue_analysis.tsv"),
            ("--matrix", dirs[5] / "rt_cluster_matrix.tsv"),
            ("--neighbor-matrix", dirs[4] / "neighbor_matrix.tsv"),
        ]
        for flag, path in optional_inputs:
            if path.exists() or args.dry_run:
                cmd.extend([flag, str(path)])

        if args.verbose:
            cmd.append("-v")

        success = run_step(cmd, "Step 7: 综合报告与类型分类", logger, args.dry_run) and success

    # ==================== Step 8: ncRNA共变验证 (可选) ====================
    if start_step <= 8 <= end_step and success and not args.skip_ncrna:
        if not fasta_dir:
            logger.warning("Step 8 需要 --fasta-dir，跳过ncRNA验证")
        else:
            # 检查工具
            tools_ok = True
            if not shutil.which('mafft'):
                logger.warning("MAFFT未安装")
                tools_ok = False
            if not shutil.which('RNAalifold'):
                logger.warning("RNAalifold (ViennaRNA)未安装")
                tools_ok = False

            if tools_ok or args.dry_run:
                # 确定输入文件
                ncrna_input = dirs[7] / "final_candidates.tsv"
                if not ncrna_input.exists() and not args.dry_run:
                    ncrna_input = dirs[3] / "retron_candidates.tsv"

                cmd = [
                    sys.executable, str(script_dir / "step08_ncrna_validation.py"),
                    "-i", str(ncrna_input),
                    "-f", str(fasta_dir),
                    "-o", str(dirs[8]),
                    "--upstream", str(upstream),
                    "--min-group-size", str(min_group_size),
                    "--threads", str(threads)
                ]

                # 添加系统发育文件用于分组
                phylo_file = dirs[3] / "retron_candidates.tsv"
                if phylo_file.exists() or args.dry_run:
                    cmd.extend(["--phylogeny", str(phylo_file)])

                if args.verbose:
                    cmd.append("-v")

                success = run_step(cmd, "Step 8: ncRNA共变验证", logger, args.dry_run) and success
            else:
                logger.warning("缺少必要工具，跳过Step 8 ncRNA共变验证")

    # 总结
    logger.info("\n" + "=" * 60)
    if success:
        logger.info("✓ 流程完成!")
    else:
        logger.error("✗ 部分步骤失败")
    logger.info(f"输出目录: {output_dir}")

    # 输出结构说明
    logger.info("\n输出目录结构:")
    for step_num, step_dir in dirs.items():
        if step_dir.exists():
            logger.info(f"  {step_dir.name}/")

    logger.info("=" * 60)

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
