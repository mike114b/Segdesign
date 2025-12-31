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
from cluster_analysis import cluster_analysis

def parse_args():
    parser = argparse.ArgumentParser(description='Protein sequence prediction', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--parse_multiple_chains_path', type=str, help='Path to parse_multiple_chains.py')
    parser.add_argument('--assign_fixed_chains_path', type=str, help='Path to assign_fixed_chains.py')
    parser.add_argument('--make_fixed_positions_dict_path', type=str, help='Path to make_fixed_positions_dict_path.py')
    parser.add_argument('--protein_mpnn_run_path', type=str, help='Path to protein_mpnn_run.py')
    parser.add_argument('--pdb_folder', type=str, help='Path of the PDB file')
    parser.add_argument('--output_folder', type=str, help='Folder for storing output files')
    parser.add_argument('--chain_list', type=str, default='A', help='Chain to be designed')
    parser.add_argument('--position_list', type=str, default=None, help='Redesigned sequence region')
    parser.add_argument('--num_seq_per_target', type=int, help='Number of generated sequences')
    parser.add_argument('--sampling_temp', type=float, default=0.1, help='Sampling temperature')
    parser.add_argument('--seed', type=int, default=37, help='Random seed')
    


    return parser.parse_args()


def parse_position_list(input_str):
    """
    解析输入字符串，支持两种格式：
    1. 纯数字格式（如 "1 2 3, 5 6 7"）
    2. 字母+数字范围格式（如 "A1-3, B5-6"）
    返回：按字母顺序排序后，仅包含数字的字符串。
    """
    # 1. 分割输入（按逗号分割，处理前后空格）
    parts = [part.strip() for part in input_str.split(',') if part.strip()]

    # 2. 分类处理每个部分
    result_parts = []
    for part in parts:
        # 判断是否为「字母+数字范围」格式
        if re.search(r'[A-Za-z]', part):
            # 提取字母（用于排序）和所有数字
            letter = re.findall(r'[A-Za-z]', part)[0].upper()
            nums = list(map(int, re.findall(r'\d+', part)))

            if len(nums) == 2 and '-' in part:
                # 情况A：范围，如 A1-3
                start_num, end_num = nums
                num_range = list(range(start_num, end_num + 1))
            else:
                # 情况B：单个数字，如 A5
                num_range = nums

            # 将数字列表转为字符串 "1 2 3"
            num_str = ' '.join(map(str, num_range))
            # 存储 (排序键, 数字字符串)
            result_parts.append((letter, num_str))

        else:
            # 纯数字格式，直接保留
            result_parts.append(('', part))

    # 3. 按字母顺序排序
    result_parts.sort(key=lambda x: x[0])

    # 4. 拼接结果，只取数字部分
    return ', '.join([num_str for _, num_str in result_parts])


def run_command(command, command_name):
    # 创建子进程，捕获标准输出和错误
    print('*'*10)
    print(f"Now starting to execute the command:\n{command}")
    print('*'*10)
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        shell=True,
        stderr=subprocess.STDOUT,  # 合并错误输出到标准输出
        bufsize=1,  # 行缓冲
        universal_newlines=True  # 返回字符串而非字节
    )
    # 实时打印输出的函数
    def print_output():
        for line in iter(process.stdout.readline, ''):
            # 移除行尾换行符后打印
            print(line, end='')
            sys.stdout.flush()  # 确保立即显示
        process.stdout.close()
    # 启动输出打印线程
    output_thread = threading.Thread(target=print_output)
    output_thread.daemon = True  # 主程序退出时自动结束线程
    output_thread.start()
    # 等待进程结束
    process.wait()
    # 检查退出状态
    if process.returncode != 0:
        raise RuntimeError(f"Command {command_name} execution failed，exit code: {process.returncode}")
    return

def extract_sequences(content):
    """从文件内容中提取序列数据"""
    sequences = []
    lines = content.strip().split('\n')

    current_seq = None
    for line in lines:
        if line.startswith('>'):
            if current_seq is not None:
                sequences.append(current_seq)
            parts = line[1:].split(', ')
            attrs = {}
            for part in parts:
                if '=' in part:
                    key, value = part.split('=', 1)
                    attrs[key.strip()] = value.strip()
            current_seq = {'header': line, 'attributes': attrs, 'sequence': []}
        elif current_seq is not None:
            current_seq['sequence'].append(line)

    if current_seq is not None:
        sequences.append(current_seq)

    return sequences


def filter_top_sequences(sequences, percent):
    """筛选global_score最高的前n*100%序列"""
    # 只处理生成序列（跳过参考序列）
    generated_seqs = [seq for seq in sequences if 'sample' in seq['attributes']]

    # 按global_score升序排序
    generated_seqs.sort(
        key=lambda x: float(x['attributes']['global_score']),
        reverse=False
    )

    # 计算需要保留的序列数量
    n = max(1, math.ceil(len(generated_seqs) * percent))
    return generated_seqs[:n]


def filtered_sequence(file_path, output_file, percent=20):
    """主处理函数"""
    with open(file_path, 'r', encoding='utf-8') as f:
        file_content = f.read()
    # 解析序列数据
    sequences = extract_sequences(file_content)

    # 筛选top序列
    top_sequences = filter_top_sequences(sequences, percent)

    # 输出结果
    if not top_sequences:
        print("No sequences found")
        return

    print(f"Top {percent}% sequences (global_score descending):{output_file}")
    with open(output_file, 'a+') as f:
        f.truncate(0)
        for seq in top_sequences:
            f.write(seq['header']+'\n')
            f.write(''.join(seq['sequence'])+'\n')
    return

#批量处理
def batch_processing(input_folder, output_folder, percent=20):
    if not os.path.exists(output_folder):  ##新建文件夹
        os.makedirs(output_folder, exist_ok=True)
    files = os.listdir(input_folder)

    for file in files:
        #filename = file.rsplit('.')[0]
        file_path = os.path.join(input_folder, file)
        out_path = os.path.join(output_folder, file)
        filtered_sequence(file_path, out_path, percent)
    return


def get_start_end(input_str):
    """
    提取输入中的开始数字和结束数字
    :param input_str: 输入字符串（格式：A1-3 / A6 / 1 2 3 等）
    :return: (start_num, end_num) 开始数字和结束数字
    """
    # 情况1：空格分隔的连续数字（如1 2 3）
    if " " in input_str:
        num_list = [int(num) for num in input_str.split()]  # 按空格分割（支持多空格）
        return num_list[0], num_list[-1]

    # 情况2：连字符分隔格式（如A1-3、1-5）
    elif "-" in input_str:
        match = re.match(r"^[A-Za-z]*(\d+)-(\d+)$", input_str)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            return start, end

    # 情况3：字母+数字（如A6）或纯数字（如6）
    else:
        match = re.match(r"^[A-Za-z]*(\d+)$", input_str)
        if match:
            num = int(match.group(1))
            return num, num

    # 无效输入返回None（可选）
    return None, None





if __name__ == "__main__":
    args = parse_args()
    out_folder = os.path.expanduser(args.output_folder)
    parse_multiple_chains_path = os.path.expanduser(args.parse_multiple_chains_path)
    assign_fixed_chains_path = os.path.expanduser(args.assign_fixed_chains_path)
    make_fixed_positions_dict_path = os.path.expanduser(args.make_fixed_positions_dict_path)
    protein_mpnn_run_path = os.path.expanduser(args.protein_mpnn_run_path)


    if not os.path.exists(out_folder):  ##新建文件夹
        os.makedirs(out_folder, exist_ok=True)
    position_list = parse_position_list(args.position_list)
    parse_multiple_chains_input = os.path.expanduser(args.pdb_folder)
    parse_multiple_chains_output = f'{out_folder}/parsed_pdbs.jsonl'
    assign_fixed_chains_input = parse_multiple_chains_output
    assign_fixed_chains_output = f'{out_folder}/assigned_pdbs.jsonl'
    make_fixed_positions_dict_input = parse_multiple_chains_output
    make_fixed_positions_dict_output = f'{out_folder}/fixed_pdbs.jsonl'

    command_parse_multiple_chains = f"""
    python {shlex.quote(parse_multiple_chains_path)} \
    --input_path={shlex.quote(parse_multiple_chains_input)} \
    --output_path={shlex.quote(parse_multiple_chains_output)} 
    """
    command_assign_fixed_chains = f"""
    python {shlex.quote(assign_fixed_chains_path)}  \
    --input_path={shlex.quote(assign_fixed_chains_input)} \
    --output_path={shlex.quote(assign_fixed_chains_output)} \
    --chain_list {shlex.quote(args.chain_list)}
    """
    command_make_fixed_dict_positions = f"""
    python {shlex.quote(make_fixed_positions_dict_path)}  \
    --input_path={shlex.quote(make_fixed_positions_dict_input)} \
    --output_path={shlex.quote(make_fixed_positions_dict_output)} \
    --chain_list {shlex.quote(args.chain_list)} \
    --position_list {shlex.quote(position_list)} \
    --specify_non_fixed
    """
    command_protein_mpnn_run = f"""
    python {shlex.quote(protein_mpnn_run_path)}  \
    --jsonl_path {shlex.quote(parse_multiple_chains_output)} \
    --chain_id_jsonl {shlex.quote(assign_fixed_chains_output)} \
    --fixed_positions_jsonl {shlex.quote(make_fixed_positions_dict_output)} \
    --out_folder {shlex.quote(out_folder)} \
    --num_seq_per_target {shlex.quote(str(args.num_seq_per_target))} \
    --sampling_temp {shlex.quote(str(args.sampling_temp))} \
    --seed {shlex.quote(str(args.seed))}
    """
    run_command(command_parse_multiple_chains, 'parse_multiple_chains')
    run_command(command_assign_fixed_chains, 'assign_fixed_chains')
    run_command(command_make_fixed_dict_positions, 'make_fixed_dict_positions')
    run_command(command_protein_mpnn_run, 'protein_mpnn_run')
     

