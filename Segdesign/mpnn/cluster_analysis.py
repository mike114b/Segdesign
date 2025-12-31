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
        cluster_rep: 代表序列FASTA文件路径
    """
    # 确保输出前缀是绝对路径
    output_prefix_abs = output_prefix.resolve()
    input_fasta_abs = input_fasta.resolve()
    
    # 确保输出目录存在
    output_prefix_abs.parent.mkdir(parents=True, exist_ok=True)
    
    # 使用骨架文件夹作为工作目录，不再使用临时目录
    import shutil
    original_cwd = os.getcwd()
    
    try:
        # 切换到骨架文件夹
        os.chdir(output_prefix_abs.parent)
        
        # 在骨架文件夹中执行mmseqs
        cmd = [
            mmseqs_path, "easy-cluster",
            str(input_fasta_abs.name),  # 使用输入文件的绝对路径
            output_prefix_abs.name,     # 使用输出前缀的名称
            "tmp_mmseqs",               # mmseqs临时文件夹
            "--threads", f"{threads}",
            "--min-seq-id", f"{min_seq_id}",
            "--cov-mode", f"{cov_mode}",
            "-c", f"{coverage}",
        ]

        print(f"运行 MMseqs2: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # 检查结果文件
        cluster_rep = output_prefix_abs.parent / f"{output_prefix_abs.name}_rep_seq.fasta"
        if not cluster_rep.exists():
            print(f"错误: 代表序列文件 {cluster_rep} 未生成", file=sys.stderr)
            sys.exit(1)
        
        # 重命名为更清晰的名称
        final_cluster_rep = output_prefix_abs.parent / "cluster_output_rep_seq.fasta"
        shutil.move(str(cluster_rep), str(final_cluster_rep))
        
        # 清理mmseqs临时文件夹
        temp_mmseqs_dir = output_prefix_abs.parent / "tmp_mmseqs"
        if temp_mmseqs_dir.exists():
            shutil.rmtree(temp_mmseqs_dir)
        
        return final_cluster_rep
        
    except Exception as e:
        print(f"MMseqs2 执行失败: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        # 恢复原始工作目录
        os.chdir(original_cwd)



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
        output_fasta: Path
):
    """
    输出代表序列的原始完整序列 FASTA（保持原始CSV文件的顺序）

    参数:
        orig_csv: 原始CSV文件
        cluster_rep: 聚类代表序列FASTA
        sub_to_orig: 子序列ID -> 原始序列ID
        output_fasta: 输出的代表序列 FASTA 文件
    """
    # 读取原始CSV文件
    df_orig = pd.read_csv(orig_csv)
    
    # 获取聚类代表序列ID集合
    rep_id_set = set(record.id for record in SeqIO.parse(cluster_rep, "fasta"))
    
    # 确保输出文件名以.fa结尾
    if not output_fasta.suffix.lower() in ['.fa', '.fasta']:
        output_fasta = output_fasta.with_suffix('.fa')
    
    # 遍历原始CSV文件的行（保持原始顺序）
    representative_count = 0
    with open(output_fasta, 'w') as f:
        for idx, row in df_orig.iterrows():
            orig_id = str(row['index'])
            sequence = str(row['sequence'])
            
            # 检查这个ID是否在代表序列中
            if orig_id in sub_to_orig:
                rep_id = sub_to_orig[orig_id]
                if rep_id in rep_id_set:
                    # 如果是代表序列，按原始顺序输出
                    f.write(f'>{orig_id}\n')
                    f.write(f'{sequence}\n')
                    representative_count += 1
    
    if representative_count > 0:
        print(f"代表序列输出完成: {representative_count} 条序列 -> {output_fasta}")
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

    # 确保input_file是Path对象
    input_file_path = Path(input_file) if isinstance(input_file, str) else input_file
    
    # 在骨架文件夹中创建子区域文件
    subregion_fasta = output_folder / "subregion_sequences.fasta"
    cluster_prefix = output_folder / "cluster_output"
    
    print(f"📁 输出到骨架文件夹: {output_folder}")
    
    sub_to_orig = extract_subregions(input_file_path, subregion_fasta, start, end)
    
    # 2. 运行聚类
    cluster_rep = run_mmseqs_cluster(
        subregion_fasta,
        cluster_prefix,
        threads=threads,
        min_seq_id=min_seq_id,
        cov_mode=cov_mode,
        coverage=coverage,
        mmseqs_path=mmseqs_path
    )

    # 3. 解析聚类结果 - 从filename参数提取骨架名称
    # filename应该是类似 "Dusp4_A_2.fa" 的格式
    skeleton_name = filename.replace('.fa', '').replace('.fasta', '')
    output_path = Path(os.path.join(work_directory, f'{skeleton_name}.fa'))
    
    # 确保输出目录存在
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # 4. 输出原始代表序列
    
    # 对于CSV输入，输出结果也需要特殊处理
    if input_file_path.suffix.lower() == '.csv':
        # 对于CSV输入，我们输出处理过的聚类结果到FASTA格式
        output_representative_sequences_from_csv(
            input_file_path, cluster_rep, sub_to_orig, output_path
        )
    else:
        # 对于FASTA输入，使用原来的处理方式
        output_representative_sequences(
            input_file_path, cluster_rep, sub_to_orig, output_path
        )
    
    print(f"聚类分析完成！")
    print(f"📊 子区域文件: {subregion_fasta}")
    print(f"📊 聚类代表序列: {cluster_rep}")
    print(f"📊 最终代表序列: {output_path}")
    return output_path



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
    # 创建result文件夹在mpnn_out目录下
    results_folder = os.path.join(output_folder, 'result')
    if not os.path.exists(results_folder):
        os.makedirs(results_folder, exist_ok=True)
    
    # 创建cluster_data文件夹在mpnn_out目录下，用于保存聚类分析相关数据
    cluster_data_folder = os.path.join(output_folder, 'cluster_data')
    if not os.path.exists(cluster_data_folder):
        os.makedirs(cluster_data_folder, exist_ok=True)
    
    folder = Path(input_folder)
    filenames = [file.name for file in folder.glob(f"*.csv") if file.is_file()]
    #filenames = os.listdir(input_folder)
    for filename in filenames:
        file_name = filename.rsplit('.')[0]
        file_path = os.path.join(input_folder, filename)

        # 确定输出文件名（保持原始骨架文件格式）
        # 从CSV文件名提取骨架ID，例如：top_90.0%_Dusp4_A_2.csv -> Dusp4_A_2
        if filename.startswith('top_') and filename.endswith('.csv'):
            # 移除 'top_90.0%_' 前缀和 '.csv' 后缀
            skeleton_name = filename[10:-4]
        elif filename.startswith('mpnn_') and filename.endswith('.csv'):
            skeleton_name = filename[5:-4]  # 移除 'mpnn_' 前缀和 '.csv' 后缀
        else:
            skeleton_name = file_name
        
        # 为当前骨架创建独立的子文件夹
        skeleton_folder = os.path.join(cluster_data_folder, skeleton_name)
        if not os.path.exists(skeleton_folder):
            os.makedirs(skeleton_folder, exist_ok=True)

        # 使用骨架文件夹作为输出目录，用于保存聚类分析的中间文件
        output_folder_path = Path(skeleton_folder)
        
        # 输出文件名
        output_filename = f"{skeleton_name}.fa"
        work_directory = results_folder  # 最终代表序列仍然输出到result目录
        
        comprehensive(
            input_file=file_path,
            output_folder=output_folder_path,  # 聚类中间文件输出到骨架特定文件夹
            filename=output_filename,  # 使用修改后的文件名
            work_directory= work_directory,  # 最终代表序列输出到result目录
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
