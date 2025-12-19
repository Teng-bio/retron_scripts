#!/bin/bash
# 构建 Retron RT 特异 HMM 模型
# 用法: bash build_retron_hmm.sh

set -e

REFERENCE_DIR="/home/teng/claude_code/retron/database/rt_reference"
OUTPUT_DIR="/home/teng/claude_code/retron/database/retron_hmm"
THREADS=8

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "Step 1: 合并参考序列"
cat "$REFERENCE_DIR"/*.fasta > reference_rt.fasta
echo "  序列数: $(grep -c '^>' reference_rt.fasta)"

echo "Step 2: MAFFT L-INS-i 比对"
mafft --localpair --maxiterate 1000 --thread $THREADS reference_rt.fasta > reference_rt_aligned.fasta

echo "Step 3: 转换为 Stockholm 格式"
python3 -c "
from Bio import SeqIO
seqs = list(SeqIO.parse('reference_rt_aligned.fasta', 'fasta'))
with open('reference_rt_aligned.sto', 'w') as f:
    f.write('# STOCKHOLM 1.0\n')
    for s in seqs: f.write(f'{s.id}  {s.seq}\n')
    f.write('//\n')
print(f'  转换完成: {len(seqs)} 条序列')
"

echo "Step 4: 构建 HMM 模型"
hmmbuild Retron_RT.hmm reference_rt_aligned.sto

echo "Step 5: 压缩索引"
hmmpress Retron_RT.hmm

echo ""
echo "完成! 输出文件:"
ls -la "$OUTPUT_DIR"
echo ""
echo "测试命令:"
echo "  hmmsearch --domtblout hits.txt Retron_RT.hmm /path/to/candidates.faa"
