import threading
import sys
import subprocess
import shlex
import argparse
import re
from pathlib import Path
import os
import pandas as pd
import shutil
import math
from Bio import SeqIO
import csv

# 导入whether_pass列添加功能
try:
    from modify_mpnn_report import add_whether_pass_column
except ImportError:
    print("警告: 无法导入modify_mpnn_report模块，whether_pass列功能将不可用")
    add_whether_pass_column = None


def parse_args():
    parser = argparse.ArgumentParser(description='Protein sequence prediction and report generation', 
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--seq_folder", type=str,
                        help="Folder containing MPNN generated fasta files")
    parser.add_argument("--output_folder", type=str,
                        help="Folder for storing output files")
    parser.add_argument("--final_report_folder", type=str, default=None,
                        help="Folder for storing final mpnn_report.csv (default: same as output_folder)")
    parser.add_argument('--top_percent', type=float, default=0.2,
                        help='Filter sequences with the lowest global_score by percentage (default: 0.2 for 20%)')
    parser.add_argument('--position_list', type=str, default=None, 
                        help='Redesigned sequence region for cluster analysis')
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="MMseqs2 number of threads (default: 8)")
    parser.add_argument("--min_seq_id", type=float, default=None,
                        help="Minimum sequence similarity (default: 0.8)")
    parser.add_argument("--cov_mode", type=int, default=0,
                        help="Coverage mode (0 = bidirectional, 1 = query, default: 0)")
    parser.add_argument("-c", "--coverage", type=float, default=0.8,
                        help="Coverage threshold (default: 0.8)")
    parser.add_argument("--mmseqs_path", type=str, default="mmseqs",
                        help="mmseqs command path (default: mmseqs)")
    parser.add_argument("-s", "--sensitivity",type=float, default=4.0,
                        help="Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [4.000]")
    parser.add_argument("--generate_report", type=bool, default=True,
                        help="Generate comprehensive MPNN report")
    parser.add_argument("--rfdiffusion_report_path", type=str, default=None,
                        help="The path to rfdiffusion_report.csv. If not entered, the default path will be used: {work_dir}/rfdiffusion_report.csv")
    
    return parser.parse_args()


def extract_sequences_from_fasta(file_path):
    """
    从FASTA文件中提取序列数据，区分初始序列和生成序列
    """
    sequences = []
    try:
        with open(file_path, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                # 解析头部信息
                header = record.description
                
                # 提取属性
                attributes = {}
                header_parts = header.split(', ')
                for part in header_parts:
                    if '=' in part:
                        key, value = part.split('=', 1)
                        attributes[key.strip()] = value.strip()
                
                sequence_data = {
                    'header': header,
                    'attributes': attributes,
                    'sequence': str(record.seq),
                    'id': record.id
                }
                sequences.append(sequence_data)
    except Exception as e:
        print(f"读取文件 {file_path} 时出错：{e}")
        return []
    
    return sequences


def natural_sort_key(filename):
    """生成自然排序的key：将文件名拆分为字符串和数字部分，数字转整数"""
    parts = re.split(r'(\d+)', os.path.splitext(filename)[0])
    key = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part)
    return key


def load_backbone_data_from_rfdiffusion(working_dir, rfdiffusion_report_path = None):
    """
    从rfdiffusion_report.csv中加载骨架数据
    """
    backbone_data = {}
    try:
        # 构建rfdiffusion_report.csv的完整路径
        if rfdiffusion_report_path == 'None' or rfdiffusion_report_path is None:
            rf_report_path = os.path.join(working_dir, 'rfdiffusion_report.csv')
        else:
            rf_report_path = rfdiffusion_report_path
        #if not os.path.exists(rf_report_path):
            # 尝试其他可能的路径
           # rf_report_path = os.path.join(working_dir, 'rfdiffusion_out', 'rfdiffusion_report.csv')
        
        if os.path.exists(rf_report_path):
            df_rf = pd.read_csv(rf_report_path)
            for _, row in df_rf.iterrows():
                backbone_index = row['index']
                backbone_data[backbone_index] = {
                    'ss8': row.get('design_ss8', ''),
                    'ss3': row.get('design_ss3', ''),
                    'H_prop': row.get('H_prop', 0.0),
                    'E_prop': row.get('E_prop', 0.0),
                    'C_prop': row.get('C_prop', 0.0),
                    'backbone': row.get('backbone', ''),
                    'success_backbone': row.get('success_backbone', ''),
                    'Success': row.get('Success', '')
                }
            print(f"已加载 {len(backbone_data)} 个骨架的数据")
        else:
            print(f"警告：找不到rfdiffusion_report.csv文件: {rf_report_path}")
            
    except Exception as e:
        print(f"读取rfdiffusion_report.csv时出错：{e}")
    
    return backbone_data


def get_design_region_positions():
    """
    获取设计区域的位置信息（346-394）
    """
    # 这里可以根据实际配置文件读取，暂时硬编码
    return 346, 394


def generate_csv_for_fasta(seq_file_path, output_folder, fa_filename, working_dir, rfdiffusion_report_path = None):
    """
    为单个FASTA文件生成CSV文件，包含完整的骨架信息和MPNN数据
    """
    print(f"处理文件：{fa_filename}")
    
    # 提取所有序列
    sequences = extract_sequences_from_fasta(seq_file_path)
    
    if not sequences:
        print(f"文件 {fa_filename} 中没有找到有效序列")
        return None
    
    # 第一个序列是初始序列，后续是生成序列
    generated_sequences = sequences[1:] if len(sequences) > 1 else []
    
    if not generated_sequences:
        print(f"文件 {fa_filename} 中没有找到生成序列")
        return None
    
    # 从rfdiffusion_report.csv加载骨架数据
    backbone_data = load_backbone_data_from_rfdiffusion(working_dir, rfdiffusion_report_path)
    #print(f"backbone data: {backbone_data}")
    
    # 提取骨架ID（从文件名如"Dusp4_A_2"）
    backbone_id = fa_filename.replace('.fa', '')
    
    # 获取设计区域位置
    design_start, design_end = get_design_region_positions()
    
    # 获取对应的骨架数据
    backbone_info = backbone_data.get(backbone_id, {
        'ss8': '',
        'ss3': '',
        'H_prop': 0.0,
        'E_prop': 0.0,
        'C_prop': 0.0,
        'backbone': '',
        'success_backbone': '',
        'Success': ''
    })
    #print(f"backbone info: {backbone_info}")
    
    # 准备CSV数据
    csv_data = []
    
    for idx, seq_data in enumerate(generated_sequences):
        # 提取MPNN属性
        score = float(seq_data['attributes'].get('score', '0.0'))
        global_score = float(seq_data['attributes'].get('global_score', '0.0'))
        
        # 计算设计区域序列（从设计区域位置提取）
        full_sequence = seq_data['sequence']
        if len(full_sequence) >= design_end:
            design_region = full_sequence[design_start-1:design_end]  # Python索引从0开始
        else:
            design_region = full_sequence
        
        csv_row = {
            'index': f"{backbone_id}_mpnn_{idx}",
            'backbone': backbone_id,
            'ss8': backbone_info['ss8'],
            'ss3': backbone_info['ss3'],
            'H_prop': backbone_info['H_prop'],
            'E_prop': backbone_info['E_prop'],
            'C_prop': backbone_info['C_prop'],
            'backbone_pdb': backbone_info['success_backbone'] if backbone_info['success_backbone'] != '-' else backbone_info['backbone'],
            'score': score,
            'global_score': global_score,
            'region': design_region,
            'sequence': full_sequence
        }
        csv_data.append(csv_row)
    
    # 生成CSV文件名
    csv_filename = f"mpnn_{backbone_id}.csv"
    csv_path = os.path.join(output_folder, csv_filename)
    
    # 保存CSV文件
    df = pd.DataFrame(csv_data)
    df.to_csv(csv_path, index=False)
    
    print(f"已生成CSV文件：{csv_filename}，包含 {len(csv_data)} 个序列")
    return csv_path, csv_data


def filter_top_sequences(csv_data, top_percent):
    """
    根据global_score筛选最低的top_percent百分比序列（保持原始index顺序）
    """
    if not csv_data:
        return []
    
    # 按global_score排序（升序，数值越低越好）以找出需要保留的序列
    sorted_data = sorted(csv_data, key=lambda x: x['global_score'])
    
    # 计算需要保留的序列数量
    total_sequences = len(sorted_data)
    n = max(1, math.ceil(total_sequences * top_percent))
    
    # 获取需要保留的序列的index（按原始顺序）
    top_indices = {seq['index'] for seq in sorted_data[:n]}
    
    # 按原始顺序返回筛选后的序列
    filtered_sequences = [seq for seq in csv_data if seq['index'] in top_indices]
    
    return filtered_sequences


def process_all_fasta_files(seq_folder, output_folder, top_percent, rfdiffusion_report_path = None):
    """
    处理所有FASTA文件并生成相应的CSV文件
    """
    print(f"开始处理FASTA文件...")
    print(f"输入文件夹：{seq_folder}")
    print(f"输出文件夹：{output_folder}")
    
    # 获取工作目录（假设seq_folder在工作目录下的mpnn_out/seqs）
    working_dir = output_folder.rsplit('/', 1)[0]
    print(f"工作目录：{working_dir}")
    
    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)
    
    # 获取所有FASTA文件
    fa_files = sorted([f for f in os.listdir(seq_folder) if f.endswith('.fa')], 
                     key=natural_sort_key)
    
    if not fa_files:
        print(f"在文件夹 {seq_folder} 中没有找到FASTA文件")
        return [], []
    
    print(f"找到 {len(fa_files)} 个FASTA文件")
    
    # 创建seqs_csv文件夹
    seqs_csv_folder = os.path.join(output_folder, 'seqs_csv')
    os.makedirs(seqs_csv_folder, exist_ok=True)
    
    # 处理每个FASTA文件
    all_csv_data = []
    generated_files = []
    
    # 创建top_percent文件夹
    top_percent_str = f"{top_percent*100:.1f}%"
    top_folder = os.path.join(output_folder, f'top_{top_percent_str}')
    os.makedirs(top_folder, exist_ok=True)
    
    # 对每个FASTA文件独立进行筛选
    top_generated_files = []
    
    for fa_file in fa_files:
        fa_file_path = os.path.join(seq_folder, fa_file)
        
        result = generate_csv_for_fasta(fa_file_path, seqs_csv_folder, fa_file, working_dir, rfdiffusion_report_path)
        if result:
            csv_path, csv_data = result
            generated_files.append(csv_path)
            all_csv_data.extend(csv_data)
            
            # 对当前文件的序列进行独立筛选
            top_sequences_current = filter_top_sequences(csv_data, top_percent)
            
            if top_sequences_current:
                base_name = os.path.splitext(fa_file)[0]
                top_csv_filename = f"top_mpnn_{base_name}.csv"
                top_csv_path = os.path.join(top_folder, top_csv_filename)
                
                df = pd.DataFrame(top_sequences_current)
                df.to_csv(top_csv_path, index=False)
                top_generated_files.append(top_csv_path)
                
                print(f"已生成Top序列CSV文件：{top_csv_filename}，包含 {len(top_sequences_current)} 个序列")
            else:
                print(f"文件 {fa_file} 中没有符合筛选条件的序列")
    
    print(f"Top序列已保存到文件夹：{top_folder}")
    
    return generated_files, top_generated_files


def get_start_end(input_str):
    """
    提取输入中的开始数字和结束数字
    """
    if " " in input_str:
        num_list = [int(num) for num in input_str.split()]
        return num_list[0], num_list[-1]
    elif "-" in input_str:
        match = re.match(r"^[A-Za-z]*(\d+)-(\d+)$", input_str)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            return start, end
    else:
        match = re.match(r"^[A-Za-z]*(\d+)$", input_str)
        if match:
            num = int(match.group(1))
            return num, num
    return None, None


def generate_final_mpnn_report(output_folder, top_percent, position_list, final_report_folder=None):
    """
    生成最终的mpnn_report.csv文件，包含所有序列
    
    参数:
        output_folder: 输出文件夹
        top_percent: top筛选百分比
        final_report_folder: 最终报告输出文件夹（默认为output_folder）
    """
    print("生成最终的MPNN报告（包含所有序列）...")
    segment = position_list
    if position_list and position_list[0].isalpha():  # 检查字符串非空且首字符是字母
        segment = position_list[1:]  # 删除首字符

    # 确定最终报告输出路径
    if final_report_folder is None:
        final_report_folder = output_folder.rsplit('/', 1)[0]
    
    # 只读取所有原始序列CSV文件
    seqs_csv_folder = os.path.join(output_folder, 'seqs_csv')
    top_percent_str = f"{top_percent*100:.1f}%"
    top_folder = os.path.join(output_folder, f'top_{top_percent_str}')
    
    # 获取所有序列的index集合（用于标记是否为Top序列）
    top_sequence_indices = set()
    if os.path.exists(top_folder):
        top_csv_files = [f for f in os.listdir(top_folder) if f.endswith('.csv')]
        for csv_file in top_csv_files:
            csv_path = os.path.join(top_folder, csv_file)
            df_top = pd.read_csv(csv_path)
            top_sequence_indices.update(df_top['index'].tolist())
    
    report_data = []
    
    # 处理所有原始序列CSV文件
    if os.path.exists(seqs_csv_folder):
        csv_files = [f for f in os.listdir(seqs_csv_folder) if f.endswith('.csv')]
        
        for csv_file in csv_files:
            csv_path = os.path.join(seqs_csv_folder, csv_file)
            df = pd.read_csv(csv_path)
            
            for _, row in df.iterrows():
                # 检查该序列是否在Top筛选中
                is_top_sequence = row['index'] in top_sequence_indices
                
                report_entry = {
                    'index': row['index'],
                    'backbone': row.get('backbone', ''),
                    'segment': segment,
                    'ss8': row.get('ss8', ''),
                    'ss3': row.get('ss3', ''),
                    'H_prop': row.get('H_prop', ''),
                    'E_prop': row.get('E_prop', ''),
                    'C_prop': row.get('C_prop', ''),
                    'backbone_pdb': row.get('backbone_pdb', ''),
                    'score': row['score'],
                    'global_core': row['global_score'],
                    'region': row.get('region', ''),
                    'sequence': row['sequence']
                }
                report_data.append(report_entry)
    
    # 生成最终报告
    if report_data:
        final_report_path = os.path.join(final_report_folder, 'mpnn_report.csv')
        df_final = pd.DataFrame(report_data)
        df_final.to_csv(final_report_path, index=False)
        
        print(f"最终MPNN报告已生成：{final_report_path}")
        print(f"包含 {len(report_data)} 条记录")
        
        # 添加whether_pass列（如果聚类分析已完成）
        if add_whether_pass_column is not None:
            try:
                result_folder = os.path.join(output_folder, 'results')
                if os.path.exists(result_folder):
                    print("🔄 开始添加whether_pass列...")
                    add_whether_pass_column(final_report_path, result_folder)
                    print("✅ whether_pass列添加成功")
                else:
                    print("ℹ️  未找到聚类结果文件夹，跳过whether_pass列添加")
            except Exception as e:
                print(f"⚠️  添加whether_pass列时出错: {e}")
        
        return final_report_path
    else:
        print("没有数据生成最终报告")
        return None


if __name__ == "__main__":
    args = parse_args()
    seq_folder = os.path.expanduser(args.seq_folder)
    output_folder = os.path.expanduser(args.output_folder)
    top_percent = args.top_percent
    rfdiffusion_report_path = args.rfdiffusion_report_path
    position_list = args.position_list
    
    print("=== MPNN序列处理和报告生成 ===")
    print(f"输入序列文件夹: {seq_folder}")
    print(f"输出文件夹: {output_folder}")
    print(f"Top筛选百分比: {top_percent*100:.1f}%")
    
    # 处理所有FASTA文件并生成CSV
    all_csv_files, top_csv_files = process_all_fasta_files(seq_folder, output_folder, top_percent, rfdiffusion_report_path)
    
    if all_csv_files:
        print(f"\n成功处理 {len(all_csv_files)} 个CSV文件")
        if top_csv_files:
            print(f"成功生成 {len(top_csv_files)} 个Top序列CSV文件")


        
        # 如果提供了args.min_seq_id，进行聚类分析
        if args.min_seq_id and top_csv_files:
            print(f"\n开始聚类分析...")
            threads = args.threads
            min_seq_id = args.min_seq_id
            cov_mode = args.cov_mode
            coverage = args.coverage
            mmseqs_path = args.mmseqs_path
            sensitivity = args.sensitivity
            
            start, end = get_start_end(position_list)
            if start is not None and end is not None:
                # 创建result文件夹在mpnn_out目录下
                results_folder = os.path.join(output_folder, 'results')
                if not os.path.exists(results_folder):
                    os.makedirs(results_folder, exist_ok=True)

                # 对每个top CSV文件进行聚类分析
                for top_csv_file in top_csv_files:
                    print(f"对文件 {os.path.basename(top_csv_file)} 进行聚类分析...")

                    # 从文件名提取骨架名称，用于匹配cluster_analysis.py中的逻辑
                    base_name = os.path.splitext(os.path.basename(top_csv_file))[0]
                    if base_name.startswith('top_mpnn_'):
                        skeleton_name = base_name.split('_',2)[-1]
                    else:
                        skeleton_name = base_name

                    # 创建骨架特定的输出文件夹
                    skeleton_folder = os.path.join(output_folder, 'cluster_data', skeleton_name)
                    os.makedirs(skeleton_folder, exist_ok=True)

                    # 直接调用cluster_analysis.py中的comprehensive函数
                    try:
                        from cluster_analysis import comprehensive
                        
                        # 调用comprehensive函数
                        comprehensive(
                            input_file=top_csv_file,
                            output_folder=Path(skeleton_folder),
                            filename=f"{skeleton_name}.fa",  # 使用骨架名称
                            work_directory=results_folder,
                            start=start,
                            end=end,
                            threads=threads,
                            min_seq_id=min_seq_id,
                            cov_mode=cov_mode,
                            coverage=coverage,
                            mmseqs_path=mmseqs_path
                        )
                        
                        print(f"聚类分析成功完成")
                        print(f"输出文件保存在: {results_folder}")
                        
                    except Exception as e:
                        print(f"聚类分析失败: {e}")
                        import traceback
                        traceback.print_exc()

                    print(f"聚类分析完成")

        # 生成最终报告
        if args.generate_report:
            final_report_path = generate_final_mpnn_report(output_folder, top_percent, position_list, args.final_report_folder)
            if final_report_path:
                print(f"\n[SUCCESS] 完整MPNN报告生成完成！")
                print(f"[OUTPUT] 主要输出文件:")
                print(f"   - 原始序列CSV: {output_folder}/seqs_csv/")
                print(f"   - Top序列CSV: {output_folder}/top_{top_percent * 100:.1f}%/")
                print(f"   - 最终报告: {final_report_path}")
    else:
        print("[ERROR] 没有成功处理任何文件")
        sys.exit(1)