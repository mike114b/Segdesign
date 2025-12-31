'''
Converts .sto alignment format to .a2m format.
'''
import argparse
import numpy as np
from Bio import SeqIO
from evcouplings.align.alignment import Alignment
from evcouplings.align.protocol import modify_alignment



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--target_seq_file', type=str, help='input filepath for the target sequence in fasta')
    parser.add_argument('--sto_alignment_file', type=str, help='input filepath for .sto')
    parser.add_argument('--output_prefix', type=str, help='output filepath prefix')
    parser.add_argument('--minimum_sequence_coverage', type=int, default=50, help='Minimum sequence coverage (percentage)')
    parser.add_argument('--minimum_column_coverage', type=int, default=70, help='Minimum_column_coverage (percentage)')

    return parser.parse_args()




def sto_to_a2m(target_seq_file, sto_alignment_file, output_prefix, minimum_sequence_coverage=50, minimum_column_coverage=70):
    def read_fasta(filename, return_ids=False):
        records = SeqIO.parse(filename, 'fasta')
        seqs = list()
        ids = list()
        for record in records:
            seqs.append(str(record.seq))
            ids.append(str(record.id))
        if return_ids:
            return seqs, ids
        else:
            return seqs

    #args = parse_args()
    with open(sto_alignment_file) as a:
        ali_raw = Alignment.from_file(a, "stockholm")

    # 聚焦目标序列区域（去除空位列）
    focus_cols = np.array([c != "-" for c in ali_raw[0]])  # 创建目标序列非空位掩码
    focus_ali = ali_raw.select(columns=focus_cols)  # 选择有效列
    target_seq, target_id = read_fasta(target_seq_file, return_ids=True)
    assert len(target_seq) == 1, 'more than 1 target seq'  # 验证单序列
    target_seq = target_seq[0]
    target_id = target_id[0]

    # 验证长度匹配
    assert len(target_seq) == len(focus_ali[0]), (
        f'{len(focus_cols)} focus cols, expected {len(target_seq)}')

    target_seq_index = 0  # 目标序列索引（第一行）
    region_start = 0  # 起始位置（完整序列）

    kwargs = {
        'prefix': output_prefix,
        'seqid_filter': None,  # 无序列ID过滤
        'hhfilter': None,  # 无HHfilter过滤
        'minimum_sequence_coverage': minimum_sequence_coverage,  # 序列最小覆盖率50%
        'minimum_column_coverage': minimum_column_coverage,  # 列最小覆盖率70%
        'compute_num_effective_seqs': False,  # 不计算有效序列数
        'theta': 0.8,  # 序列聚类阈值（未使用）
    }
    #print('focus_ali:', focus_ali)
    #print("=" * 50)
    #print("🔍 调试：Alignment 对象（focus_ali）详细信息")
    #print(f"  1. 对齐序列总数：{len(focus_ali)}")  # 关键：是否有有效序列（需 ≥1）
    #print(f"  2. 所有序列ID列表：{focus_ali.ids}")  # 关键：目标序列ID是否在其中
    #print(f"  3. 序列长度（对齐后）：{focus_ali.length}")  # 对齐后的列数（需 ≥1）
    #print(f"  4. 目标序列ID（传入的 target_id）：{target_id}")
    #print(f"  5. 硬编码的目标索引（target_seq_index）：{target_seq_index}")
    #print("=" * 50)
    mod_outcfg, ali = modify_alignment(
        focus_ali, target_seq_index, target_id, region_start, **kwargs
    )

    return

def main():
    args = parse_args()
    sto_to_a2m(args.target_seq_file, args.sto_alignment_file, args.output_prefix)
    return

if __name__ == "__main__":
    main()
