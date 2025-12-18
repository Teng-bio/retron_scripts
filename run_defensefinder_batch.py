#!/usr/bin/env python3
"""
批量运行DefenseFinder的独立脚本

此脚本用于在DefenseFinder专属conda环境中运行。
它会遍历step06生成的所有基因簇，对每个基因簇的邻近蛋白运行DefenseFinder。

使用方法:
  # 激活DefenseFinder环境
  conda activate defensefinder

  # 运行批量分析
  python run_defensefinder_batch.py -i 06_neighbors/clusters -o defensefinder_results

  # 然后在主环境中运行step07整合结果
  conda activate main_env
  python step07_defense_island.py -i 06_neighbors --df-results defensefinder_results -o 07_defense_island
"""

import argparse
import subprocess
import sys
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


def run_defensefinder(fasta_file, output_dir):
    """对单个fasta文件运行DefenseFinder"""
    try:
        cmd = [
            'defense-finder', 'run',
            str(fasta_file),
            '-o', str(output_dir),
            '--db-type', 'unordered'
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300
        )

        # 检查是否生成了结果文件
        genes_file = output_dir / "defense_finder_genes.tsv"
        return genes_file.exists()

    except subprocess.TimeoutExpired:
        return False
    except Exception as e:
        print(f"  错误: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="批量运行DefenseFinder (在专属环境中使用)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  conda activate defensefinder
  python run_defensefinder_batch.py -i 06_neighbors/clusters -o defensefinder_results

注意:
  - 运行前请先执行: defense-finder update
  - 结果将保存在输出目录的子目录中，每个基因簇一个目录
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='基因簇目录 (06_neighbors/clusters)')
    parser.add_argument('-o', '--output', default='defensefinder_results',
                        help='输出目录')
    parser.add_argument('-j', '--jobs', type=int, default=1,
                        help='并行任务数 (默认: 1)')

    args = parser.parse_args()

    # 检查DefenseFinder
    if not shutil.which('defense-finder'):
        print("错误: DefenseFinder未安装或不在PATH中", file=sys.stderr)
        print("请确保已激活正确的conda环境", file=sys.stderr)
        sys.exit(1)

    clusters_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 查找所有基因簇
    cluster_dirs = [d for d in clusters_dir.iterdir() if d.is_dir()]
    print(f"找到 {len(cluster_dirs)} 个基因簇")

    # 过滤出有蛋白文件的基因簇
    tasks = []
    for cluster_dir in cluster_dirs:
        fasta_file = cluster_dir / "neighbor_proteins.fasta"
        if fasta_file.exists():
            cluster_id = cluster_dir.name
            cluster_output = output_dir / cluster_id
            cluster_output.mkdir(parents=True, exist_ok=True)
            tasks.append((cluster_id, fasta_file, cluster_output))

    print(f"有效基因簇: {len(tasks)}")

    # 运行DefenseFinder
    success_count = 0
    defense_found = 0

    if args.jobs > 1:
        # 并行运行
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(run_defensefinder, fasta, out): cid
                for cid, fasta, out in tasks
            }

            for i, future in enumerate(as_completed(futures)):
                cluster_id = futures[future]
                try:
                    has_defense = future.result()
                    success_count += 1
                    if has_defense:
                        defense_found += 1
                except Exception as e:
                    print(f"  {cluster_id}: 失败 - {e}")

                if (i + 1) % 50 == 0:
                    print(f"进度: {i + 1}/{len(tasks)}")
    else:
        # 串行运行
        for i, (cluster_id, fasta_file, cluster_output) in enumerate(tasks):
            has_defense = run_defensefinder(fasta_file, cluster_output)
            success_count += 1
            if has_defense:
                defense_found += 1
                print(f"  [{i+1}/{len(tasks)}] {cluster_id}: 发现防御系统")

            if (i + 1) % 50 == 0:
                print(f"进度: {i + 1}/{len(tasks)}")

    # 总结
    print("\n" + "=" * 50)
    print("DefenseFinder批量分析完成!")
    print(f"  处理基因簇: {success_count}/{len(tasks)}")
    print(f"  发现防御系统: {defense_found}")
    print(f"  输出目录: {output_dir}")
    print("=" * 50)
    print("\n下一步: 在主环境中运行step07整合结果:")
    print(f"  python step07_defense_island.py -i {clusters_dir.parent} --df-results {output_dir} -o 07_defense_island")

    return 0


if __name__ == "__main__":
    sys.exit(main())
