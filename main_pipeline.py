#!/usr/bin/env python3
"""
Retron系统挖掘流程 v2.0 - 主控脚本

基于Millman et al., 2020 Cell的方法论，整合多维度筛选策略：
1. RT搜索 (DIAMOND/BLAST)
2. Motif过滤 (NAXXH + VTG)
3. RT分类与排除 (Group II/DGR/CRISPR)
4. 上游序列提取
5. RNA结构筛选 (滑动窗口 + 分支G)
6. 邻近基因提取
7. 防御岛分析 (DefenseFinder)
8. 效应蛋白检测
9. 综合报告与置信度分级
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
        description="Retron系统挖掘流程 v2.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 使用配置文件
  python main_pipeline.py -c config.yaml

  # 命令行参数
  python main_pipeline.py \\
    --search-results /path/to/complete_results.tsv \\
    --faa-dir /path/to/faa \\
    --fasta-dir /path/to/fasta \\
    --antismash-dir /path/to/antismash \\
    --output ./retron_v2_results

  # 预览命令
  python main_pipeline.py -c config.yaml --dry-run
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
    parser.add_argument('--start-step', type=int, default=1, help='起始步骤 (1-9)')
    parser.add_argument('--end-step', type=int, default=9, help='结束步骤 (1-9)')
    parser.add_argument('--skip-defense', action='store_true',
                        help='跳过防御岛分析 (Step 7)')
    parser.add_argument('--df-results',
                        help='已有的DefenseFinder结果目录 (用于Step 7)')

    # 参数
    parser.add_argument('--min-identity', type=float, default=25)
    parser.add_argument('--min-coverage', type=float, default=40)
    parser.add_argument('--upstream', type=int, default=600)
    parser.add_argument('--window', type=int, default=150)
    parser.add_argument('--min-mfe', type=float, default=-15.0)
    parser.add_argument('--require-branch-g', action='store_true')
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
    # 输出目录：命令行参数 > 配置文件 > 默认值
    output_dir = Path(args.output) if args.output else Path(get_config('paths.output_dir', './retron_v2_results'))

    # 显示关键配置
    if config:
        print(f"   配置的输出目录: {get_config('paths.output_dir', '未指定')}")
        print(f"   实际输出目录: {output_dir.absolute()}")

    # 参数配置
    min_identity = args.min_identity or get_config('search.min_identity', 25)
    min_coverage = args.min_coverage or get_config('search.min_coverage', 40)
    upstream = args.upstream or get_config('flanking.upstream', 600)
    window = args.window or get_config('rna_filter.window_size', 150)
    min_mfe = args.min_mfe or get_config('rna_filter.min_mfe', -15.0)
    require_branch_g = args.require_branch_g or get_config('rna_filter.require_branch_g', False)
    threads = args.threads or get_config('general.threads', 8)

    start_step = args.start_step or get_config('pipeline.start_step', 1)
    end_step = args.end_step or get_config('pipeline.end_step', 9)

    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging(output_dir, args.verbose)

    # 脚本目录
    script_dir = Path(__file__).parent

    # 步骤输出目录
    dirs = {
        1: output_dir / "01_search",
        2: output_dir / "02_motif_filtered",
        3: output_dir / "03_classified",
        4: output_dir / "04_flanking",
        5: output_dir / "05_rna_filtered",
        6: output_dir / "06_neighbors",
        7: output_dir / "07_defense_island",
        8: output_dir / "08_effector",
        9: output_dir / "09_report"
    }

    logger.info("=" * 60)
    logger.info("Retron系统挖掘流程 v2.0")
    logger.info("=" * 60)
    logger.info(f"输出目录: {output_dir}")
    logger.info(f"步骤范围: {start_step} - {end_step}")

    success = True

    # 根据start_step确定正确的初始输入文件
    step_outputs = {
        1: dirs[1] / "reports" / "complete_results.tsv",
        2: dirs[2] / "motif_filtered.tsv",
        3: dirs[3] / "retron_candidates.tsv",
        4: dirs[4] / "flanking_sequences.fasta",
        5: dirs[5] / "filter_results.tsv",
        6: dirs[6] / "neighbors_info.tsv",
        7: dirs[7] / "defense_analysis.tsv",
        8: dirs[8] / "effector_analysis.tsv",
    }

    # 如果从中间步骤开始，使用前一步骤的输出
    if start_step > 1:
        prev_step = start_step - 1
        if prev_step in step_outputs and step_outputs[prev_step].exists():
            current_input = step_outputs[prev_step]
            logger.info(f"从Step {start_step}开始，使用输入: {current_input}")
        else:
            current_input = search_results if search_results else step_outputs.get(1)
    else:
        current_input = search_results

    # Step 1: RT搜索
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
        # Step 1在范围内但跳过（提供了search_results）
        current_input = Path(search_results) if search_results else dirs[1] / "reports" / "complete_results.tsv"
    # 如果start_step > 1，current_input已经在前面设置好了

    # Step 2: Motif过滤
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

    # Step 3: RT分类
    if start_step <= 3 <= end_step and success:
        cmd = [
            sys.executable, str(script_dir / "step03_classify_rt.py"),
            "-i", str(current_input),
            "-o", str(dirs[3]),
            "--keep-excluded"
        ]
        success = run_step(cmd, "Step 3: RT分类与排除", logger, args.dry_run) and success
        current_input = dirs[3] / "retron_candidates.tsv"

    # Step 4: 上游序列提取
    if start_step <= 4 <= end_step and success:
        if not fasta_dir:
            logger.error("Step 4 需要 --fasta-dir")
            sys.exit(1)

        cmd = [
            sys.executable, str(script_dir / "step04_extract_flanking.py"),
            "-i", str(current_input),
            "-f", str(fasta_dir),
            "-o", str(dirs[4]),
            "--upstream", str(upstream)
        ]
        success = run_step(cmd, "Step 4: 上游序列提取", logger, args.dry_run) and success

    # Step 5: RNA结构筛选
    if start_step <= 5 <= end_step and success:
        cmd = [
            sys.executable, str(script_dir / "step05_rna_filter.py"),
            "-i", str(dirs[4] / "flanking_sequences.fasta"),
            "-o", str(dirs[5]),
            "--window", str(window),
            "--min-mfe", str(min_mfe)
        ]
        if require_branch_g:
            cmd.append("--require-branch-g")

        success = run_step(cmd, "Step 5: RNA结构筛选", logger, args.dry_run) and success
        current_input = dirs[5] / "filter_results.tsv"

    # Step 6: 邻近基因提取
    if start_step <= 6 <= end_step and success:
        if not antismash_dir:
            logger.warning("Step 6 需要 --antismash-dir，跳过")
        else:
            cmd = [
                sys.executable, str(script_dir / "step06_extract_neighbors.py"),
                "-i", str(current_input),
                "-a", str(antismash_dir),
                "-o", str(dirs[6])
            ]
            success = run_step(cmd, "Step 6: 邻近基因提取", logger, args.dry_run) and success

    # Step 7: 防御岛分析 (基于step06的基因簇)
    if start_step <= 7 <= end_step and success and not args.skip_defense:
        # 使用step06的输出目录作为输入
        step06_dir = dirs[6]
        cmd = [
            sys.executable, str(script_dir / "step07_defense_island.py"),
            "-i", str(step06_dir),
            "-o", str(dirs[7])
        ]
        # 如果提供了已有的DefenseFinder结果，使用它
        if args.df_results:
            cmd.extend(["--df-results", str(args.df_results)])
        elif not shutil.which('defense-finder'):
            logger.warning("DefenseFinder未安装，请使用 --df-results 提供已有结果，或跳过Step 7")
            logger.warning("提示: 先在DefenseFinder环境中运行 run_defensefinder_batch.py")
        success = run_step(cmd, "Step 7: 防御岛分析", logger, args.dry_run) and success

    # Step 8: 效应蛋白检测
    if start_step <= 8 <= end_step and success:
        neighbors_file = dirs[6] / "neighbors_info.tsv"
        if neighbors_file.exists() or args.dry_run:
            cmd = [
                sys.executable, str(script_dir / "step08_effector_detect.py"),
                "-i", str(current_input),
                "-n", str(neighbors_file),
                "-o", str(dirs[8])
            ]
            success = run_step(cmd, "Step 8: 效应蛋白检测", logger, args.dry_run) and success
        else:
            logger.warning("邻近基因文件不存在，跳过Step 8")

    # Step 9: 综合报告
    if start_step <= 9 <= end_step and success:
        cmd = [
            sys.executable, str(script_dir / "step09_final_report.py"),
            "-o", str(dirs[9]),
            "--search", str(dirs[1] / "reports" / "complete_results.tsv") if (dirs[1] / "reports" / "complete_results.tsv").exists() else str(search_results or ""),
            "--motif", str(dirs[2] / "motif_filtered.tsv"),
            "--classified", str(dirs[3] / "retron_candidates.tsv"),
            "--rna", str(dirs[5] / "filter_results.tsv"),
            "--neighbors", str(dirs[6] / "neighbors_info.tsv"),
            "--defense", str(dirs[7] / "defense_analysis.tsv"),
            "--effector", str(dirs[8] / "effector_analysis.tsv")
        ]
        # 只传入存在的文件
        success = run_step(cmd, "Step 9: 综合报告", logger, args.dry_run) and success

    # 总结
    logger.info("\n" + "=" * 60)
    if success:
        logger.info("✓ 流程完成!")
    else:
        logger.error("✗ 部分步骤失败")
    logger.info(f"输出目录: {output_dir}")
    logger.info("=" * 60)

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
