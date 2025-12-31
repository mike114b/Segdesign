#!/usr/bin/env python3
"""
MMseqs2 特定区域聚类工具
功能：提取序列指定区域 → 聚类 → 输出原始代表序列 FASTA
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

def arg_parser():
    parser = argparse.ArgumentParser(
        description="Perform MMseqs2 clustering on specific regions of sequences and output the original complete representative sequences"
    )
    parser.add_argument("-i", "--input_folder", required=True, type=Path,
                        help="Input FASTA folder or CSV folder")
    parser.add_argument("-o", "--output_folder", required=True, type=Path,
                        help="Output folder")
    parser.add_argument("-s", "--start", required=True, type=int,
                        help="Start position (1-based, inclusive)")
    parser.add_argument("-e", "--end", required=True, type=int,
                        help="End position (1-based, inclusive)")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of threads for MMseqs2 (default: 8)")
    parser.add_argument("--min_seq_id", type=float, default=0.5,
                        help="Minimum sequence identity (default: 0.5)")
    parser.add_argument("--cov_mode", type=int, default=0,
                        help="Coverage mode (0=bidirectional, 1=query, default: 0)")
    parser.add_argument("-c", "--coverage", type=float, default=0.8,
                        help="Coverage threshold (default: 0.8)")
    parser.add_argument("--mmseqs_path", type=str, default="mmseqs",
                        help="Path to mmseqs command (default: mmseqs)")
    return parser.parse_args()




def extract_subregions(
        input_file: Path,
        output_fasta: Path,
        start_pos: int,
        end_pos: int,
) -> Dict[str, str]:
    """
    从 FASTA 文件或 CSV 文件中提取特定区域，并记录 ID 映射关系

    参数:
        input_file: 输入 FASTA 文件或 CSV 文件
        output_fasta: 输出的子区域 FASTA 文件
        start_pos: 起始位置 (1-based, 包含)
        end_pos: 结束位置 (1-based, 包含)

    返回:
        sub_to_orig: 子序列ID -> 原始序列ID 的字典
    """
    sub_to_orig = {}
    sub_records = []
    ndx = 0
    
    # 检测文件类型
    if input_file.suffix.lower() == '.csv':
        # 处理CSV文件
        print(f"检测到CSV文件: {input_file}")
        try:
            df = pd.read_csv(input_file)
            if 'sequence' not in df.columns:
                print(f"错误: CSV文件 {input_file} 中没有 'sequence' 列", file=sys.stderr)
                sys.exit(1)
            if 'index' not in df.columns:
                print(f"错误: CSV文件 {input_file} 中没有 'index' 列", file=sys.stderr)
                sys.exit(1)
            
            # 保持CSV文件的原始顺序
            for idx, row in df.iterrows():
                orig_id = str(row['index'])
                sequence = str(row['sequence'])
                
                # 保持原始ID作为子序列ID，不使用简单数字
                sub_id = orig_id  # 直接使用原始ID
                sub_to_orig[sub_id] = orig_id

                # 提取子序列 (转换为0-based索引)
                start_idx = max(0, start_pos - 1)
                end_idx = min(len(sequence), end_pos)

                if start_idx >= end_idx:
                    print(f"警告: 序列 {orig_id} 长度 {len(sequence)} 小于指定区域，跳过", file=sys.stderr)
                    continue

                sub_seq = Seq(sequence[start_idx:end_idx])

                # 创建新记录
                sub_record = SeqRecord(
                    seq=sub_seq,
                    id=sub_id,
                    description=f""
                )
                sub_records.append(sub_record)
                
        except Exception as e:
            print(f"错误: 读取CSV文件失败: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # 处理FASTA文件
        print(f"检测到FASTA文件: {input_file}")
        for record in SeqIO.parse(input_file, "fasta"):
            orig_id = record.description
            # 创建子序列ID
            ndx += 1
            sub_id = f"{ndx}"
            sub_to_orig[sub_id] = orig_id

            # 提取子序列 (转换为0-based索引)
            start_idx = max(0, start_pos - 1)
            end_idx = min(len(record.seq), end_pos)

            if start_idx >= end_idx:
                print(f"警告: 序列 {orig_id} 长度 {len(record.seq)} 小于指定区域，跳过", file=sys.stderr)
                continue

            sub_seq = Seq(str(record.seq[start_idx:end_idx]))

            # 创建新记录
            sub_record = SeqRecord(
                seq=sub_seq,
                id=sub_id,
                description=f""
            )
            sub_records.append(sub_record)

    # 写入子序列FASTA
    with open(output_fasta, 'w') as f:
        SeqIO.write(sub_records, f, 'fasta')

    print(f"提取完成: {len(sub_records)} 条序列 -> {output_fasta}")
    return sub_to_orig


def run_mmseqs_cluster(
        input_fasta: Path,
        output_prefix: Path,
        threads: int = 8,
        min_seq_id: float = 0.5,
        cov_mode: int = 0,
        coverage: float = 0.8,
        mmseqs_path: str = "mmseqs"
) -> Path:
    """
    运行 MMseqs2 聚类

    参数:
        input_fasta: 输入FASTA文件
        output_prefix: 输出文件前缀
        threads: 线程数
        min_seq_id: 最小序列相似度
        cov_mode: 覆盖度模式
        coverage: 覆盖度阈值
        mmseqs_path: mmseqs 命令路径

    返回:
        cluster_tsv: 聚类结果 TSV 文件路径
    """
    cmd = [
        mmseqs_path, "easy-cluster",
        str(input_fasta),
        str(output_prefix),
        "tmp",
        "--threads", f"{threads}",
        "--min-seq-id", f"{min_seq_id}",
        "--cov-mode", f"{cov_mode}",
        "-c", f"{coverage}",
    ]

    print(f"运行 MMseqs2: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"MMseqs2 执行失败: {e.stderr}", file=sys.stderr)
        sys.exit(1)

    # 聚类结果文件
    cluster_rep = Path(f"{output_prefix}_rep_seq.fasta")
    if not cluster_rep.exists():
        print(f"错误: 代表序列文件 {cluster_rep} 未生成", file=sys.stderr)
        sys.exit(1)

    return cluster_rep



def output_representative_sequences(
        orig_fasta: Path,
        cluster_rep: Path,
        sub_to_orig: Dict[str, str],
        output_fasta: Path
):
    """
    输出代表序列的原始完整序列 FASTA

    参数:
        orig_fasta: 原始完整序列 FASTA
        rep_seq_map: 簇ID -> 子序列代表ID
        sub_to_orig: 子序列ID -> 原始序列ID
        output_fasta: 输出的代表序列 FASTA 文件
    """
    result_records = []
    # 加载原始序列到字典
    orig_records = {record.description: record.seq for record in SeqIO.parse(orig_fasta, "fasta")}
    #print('orig_records:', orig_records)
    rep_id_l = [record.id for record in SeqIO.parse(cluster_rep, "fasta")]
    with open(output_fasta, 'w') as f:
        f.truncate(0)
        for rep_id in rep_id_l:
            result_id = sub_to_orig[rep_id]
            # print('result_id:', result_id)
            result_seq = orig_records[result_id]
            # print(f'result_seq:', result_seq)

            # result_record = SeqRecord(
            # seq=result_seq,
            # id=result_id,
            # description=f""
            # )
            # result_records.append(result_record)
            f.write('>'+str(result_id) + '\n')
            f.write(str(result_seq) + '\n')


    #with open(output_fasta, 'w') as f:
        #SeqIO.write(result_records, f, 'fasta')


    print(f"代表序列输出完成: {len(rep_id_l)} 条序列 -> {output_fasta}")


def output_representative_sequences_from_csv(
        orig_csv: Path,
        cluster_rep: Path,
        sub_to_orig: Dict[str, str],
        output_csv: Path
):
    """
    输出代表序列的原始完整序列 CSV

    参数:
        orig_csv: 原始CSV文件
        cluster_rep: 聚类代表序列FASTA
        sub_to_orig: 子序列ID -> 原始序列ID
        output_csv: 输出的代表序列 CSV 文件
    """
    # 读取原始CSV文件
    df_orig = pd.read_csv(orig_csv)
    
    # 获取聚类代表序列ID
    rep_id_l = [record.id for record in SeqIO.parse(cluster_rep, "fasta")]
    
    # 创建代表序列的DataFrame
    representative_data = []
    for rep_id in rep_id_l:
        orig_id = sub_to_orig[rep_id]
        # 在原始CSV中找到对应的序列
        matching_rows = df_orig[df_orig['index'] == orig_id]
        if not matching_rows.empty:
            # 如果找到匹配的序列，添加到结果中
            row = matching_rows.iloc[0].copy()
            representative_data.append(row)
        else:
            print(f"警告: 在CSV中找不到序列ID {orig_id}", file=sys.stderr)
    
    # 保存结果到CSV
    if representative_data:
        df_result = pd.DataFrame(representative_data)
        # 确保输出文件名以.csv结尾
        if not output_csv.suffix.lower() == '.csv':
            output_csv = output_csv.with_suffix('.csv')
        df_result.to_csv(output_csv, index=False)
        print(f"代表序列输出完成: {len(representative_data)} 条序列 -> {output_csv}")
    else:
        print("错误: 没有找到代表序列", file=sys.stderr)

#整合
def comprehensive(
        input_file,
        output_folder,
        filename,
        work_directory,
        start: int,
        end: int,
        threads: int = 8,
        min_seq_id = 0.8,
        cov_mode = 0,
        coverage = 0.8,
        mmseqs_path = 'mmseqs'
):
    """
    将上述的子程序整合，实现输入FASTA或CSV文件，直接输出聚类结果

    参数:
        input_file: 输入 FASTA 文件或 CSV 文件
        output_folder: 输出目录
        start: 起始位置
        end: 结束位置
        threads: 线程数 (默认: 8)
        min_seq_id：最小序列相似度 (默认: 0.8)
        cov_mode: 覆盖度模式 (0=双向, 1=查询, 默认: 0)
        coverage: 覆盖度阈值 (默认: 0.8)
        mmseqs_path: mmseqs 命令路径 (默认: mmseqs)

    返回:
        聚类结果（存放在输出目录中）
    """

    tmp_sub_fasta = Path(os.path.join(output_folder,"tmp_subregion.fasta"))
    #print('tmp_sub_fasta:', tmp_sub_fasta)

    sub_to_orig = extract_subregions(input_file, tmp_sub_fasta, start, end)
    #print('sub_to_orig:', sub_to_orig)

    # 2. 运行聚类
    cluster_prefix = Path(os.path.join(output_folder,"cluster_output"))
    cluster_rep = run_mmseqs_cluster(
        tmp_sub_fasta,
        cluster_prefix,
        threads=threads,
        min_seq_id=min_seq_id,
        cov_mode=cov_mode,
        coverage=coverage,
        mmseqs_path=mmseqs_path
    )

    # 3. 解析聚类结果
    output_path = Path(os.path.join(work_directory, f'{filename}'))
    # 4. 输出原始代表序列
    
    # 对于CSV输入，输出结果也需要特殊处理
    if input_file.suffix.lower() == '.csv':
        # 对于CSV输入，我们输出处理过的聚类结果到CSV格式
        output_representative_sequences_from_csv(
            input_file, cluster_rep, sub_to_orig, output_path
        )
    else:
        # 对于FASTA输入，使用原来的处理方式
        output_representative_sequences(
            input_file, cluster_rep, sub_to_orig, output_path
        )

    # 5. 清理临时文件（可选）
    print("清理临时文件...\n")
    subprocess.run(["rm", "-rf", f"{output_folder}/tmp_subregion.fasta", f"tmp", f"{output_folder}/{cluster_prefix}*"], check=False)
    if os.path.exists(f"embedding_lookup_avx2.cc"):
        subprocess.run(["rm", "-rf", f"embedding_lookup_avx2.cc"], check=False)

    #print("\n✅ 所有步骤完成！")
    return



def cluster_analysis(
        input_folder,
        output_folder,
        start,
        end,
        threads=8,
        min_seq_id=0.8,
        cov_mode=0,
        coverage = 0.8,
        mmseqs_path = 'mmseqs'
):
    cluster_folder = os.path.join(output_folder, f'cluster_data')
    if not os.path.exists(cluster_folder):
        os.makedirs(cluster_folder, exist_ok=True)

    # 创建result文件夹在mpnn_out目录下
    results_folder = os.path.join(output_folder, 'result')
    if not os.path.exists(results_folder):
        os.makedirs(results_folder, exist_ok=True)

    folder = Path(input_folder)
    filenames = [file.name for file in folder.glob(f"*.csv") if file.is_file()]
    #filenames = os.listdir(input_folder)
    for filename in filenames:
        file_name = filename.rsplit('.')[0]
        file_path = os.path.join(input_folder, filename)
        output_folder_path = Path(os.path.join(cluster_folder, f'{file_name}'))
        
        if not os.path.exists(output_folder_path):
            os.makedirs(output_folder_path, exist_ok=True)

        # 确定输出文件名（保持原始骨架文件格式）
        if '_' in filename:
            # 处理格式：top_mpnn_Dusp4_A_2.csv -> Dusp4_A_2.fa
            parts = filename.split('_')
            if len(parts) >= 5:  # 确保有足够的部分
                base_name = f"{parts[2]}_{parts[3]}_{parts[4].split('.')[0]}"
            else:
                base_name = file_name
        else:
            base_name = file_name
        output_filename = f"{base_name}.fa"
        work_directory = results_folder
        
        comprehensive(
            input_file=file_path,
            output_folder=output_folder_path,
            filename=output_filename,  # 使用修改后的文件名
            work_directory= work_directory,
            start=start,
            end=end,
            threads = threads,
            min_seq_id = min_seq_id,
            cov_mode = cov_mode,
            coverage = coverage,
            mmseqs_path = mmseqs_path
        )
    print("\n✅ 所有步骤完成！")
    return







def main():
    args = arg_parser()
    input_folder = os.path.expanduser(args.input_folder)
    output_folder = os.path.expanduser(args.output_folder)
    start = args.start
    end = args.end
    threads = args.threads
    print('threads:', threads)
    min_seq_id = args.min_seq_id
    cov_mode = args.cov_mode
    coverage = args.coverage
    mmseqs_path = args.mmseqs_path

    cluster_analysis(
        input_folder=input_folder,
        output_folder=output_folder,
        start=start,
        end=end,
        threads=threads,
        min_seq_id=min_seq_id,
        cov_mode=cov_mode,
        coverage=coverage,
        mmseqs_path=mmseqs_path
    )
    return


if __name__ == "__main__":
    main()
