import csv
import threading
import sys
import subprocess
import shlex
import argparse
import re
from pathlib import Path
import os
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))  # rfdiffusion/目录
PARENT_DIR = os.path.dirname(ROOT_DIR)  # 上级目录（包含dssp/和rfdiffusion/）
sys.path.append(PARENT_DIR)  # 把上级目录加入搜索路径
from pdbrepair import fix_pdb_file
from dssp.dssp import run_dssp
from dssp.dsspcsv import dssp_to_csv
import pandas as pd
import shutil


def parse_args():
    parser = argparse.ArgumentParser(description='Protein backbone remodeling', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--run_inference_path', type=str, help='Path to run_inference.py')
    return parser.parse_known_args()


def run_command(command):
    # 创建子进程，捕获标准输出和错误
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
        raise RuntimeError(f"HMMER execution failed，exit code: {process.returncode}")
    return

#外部输入参数
def parse_parameter_file(filepath):
    """
    Parse parameter.txt and extract script path & parameters.
    Returns: (script_path, params_dict)
    """
    params = {}
    script_path = None

    with open(filepath, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            # Remove leading/trailing whitespace
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            # Split on first '=' only
            if '=' not in line:
                print(f"Warning: Line {line_num} skipped (no '=' found): {line}")
                continue

            key, value = line.split('=', 1)
            key = key.strip()
            value = value.strip()

            # Remove inline comments (preserve # in paths if quoted)
            if '#' in value and not (value.startswith('"') or value.startswith("'")):
                value = value.split('#')[0].strip()

            # Skip the rfdiffusion flag line
            if key.lower() == 'rfdiffusion':
                continue

            # Extract script path
            if key == 'run_inference_path':
                script_path = value
                continue

            # Store parameter
            params[key] = value

    if not script_path:
        raise ValueError("Required 'run_inference_path' not found in parameter file")

    return script_path, params

#运行rfdiffusion，生成蛋白质骨架pdb文件
def run_rfdiffusion(args, unknown):
    script_path = os.path.expanduser(args.run_inference_path)
    params = {}
    for i in range(0, len(unknown), 2):
        if unknown[i].startswith('--'):
            key = unknown[i][2:]
            if i + 1 < len(unknown):
                v = unknown[i + 1].strip()

                params[key] = v
    cmd = ['python', script_path]
    for key, value in params.items():
        print('value:',value)
        cmd.append(f"{key}={value}")

    # Print command for verification
    print("=" * 60)
    print("Executing RFdiffusion with command:")
    print(' '.join(shlex.quote(arg) for arg in cmd))
    print("=" * 60)
    print()
    # Run command
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=False,  # Show real-time output
            text=True
        )
        print("\n RFdiffusion completed successfully!")

    except subprocess.CalledProcessError as e:
        print(f"\n RFdiffusion failed with return code {e.returncode}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("\n Error: 'python' command or script not found. Check your environment.", file=sys.stderr)
        sys.exit(1)
    
    return



#修正蛋白质骨架pdb文件，使之可以被dssp识别
def fix_pdb(path,output_folder):
    folder_path, file_prefix = path.rsplit('/', 1)
    pattern = re.compile(f"^{file_prefix}_\d+\.pdb$")
    target_folder = Path(folder_path)
    matched_files = [
        file for file in target_folder.iterdir()  # 遍历当前文件夹所有条目
        if file.is_file()  # 仅保留文件（过滤文件夹）
           and pattern.match(file.name)  # 正则匹配文件名
    ]
    if not os.path.exists(output_folder):  ##新建文件夹
        os.makedirs(output_folder, exist_ok=True)
    for file_path in matched_files:
        file_name = str(file_path).rsplit('/', 1)[1]
        fix_pdb_file(file_path, f'{output_folder}/{file_name}')

    return

# 批量生成rfdiffusion创造的蛋白质骨架文件的二级结构文件（dssp）
def dssp_generation(input_folder,output_folder):
    filenames = os.listdir(input_folder)
    if not os.path.exists(output_folder):  ##新建文件夹
        os.makedirs(output_folder, exist_ok=True)
    for filename in filenames:
        file_name = filename.rsplit('.', 1)[0]
        run_dssp(f'{input_folder}/{filename}', f'{output_folder}/{file_name}.dssp')
        print(f'{file_name}.dssp generated successfully!')
    return

#从dssp文件中提取二级结构信息，存入csv文件
def dsspcsv(input_folder,output_folder):
    if not os.path.exists(output_folder):  ##新建文件夹
        os.makedirs(output_folder, exist_ok=True)
    filenames = os.listdir(input_folder)
    for filename in filenames:
        file_name = filename.rsplit('.', 1)[0]
        dssp_to_csv(f'{input_folder}/{filename}', f'{output_folder}/{file_name}.csv')
    return

#
def dssp_analyse(path, start_res, end_res, target_ss, threshold = 0.6):
    path_, name_prefix= path.rsplit('/', 1)
    print('path_:', path_)
    fix_pdb_folder = 'fix_pdb'
    dssp_folder = 'dssp_file'
    dssp_csv_folder = 'dssp_csv_file'
    print(os.path.exists(f'{path_}/{fix_pdb_folder}'))
    if not os.path.exists(f'{path_}/{fix_pdb_folder}'):
        os.makedirs(f'{path_}/{fix_pdb_folder}', exist_ok=True)
        print(f'{path_}/{fix_pdb_folder}')
    if not os.path.exists(f'{path_}/{dssp_folder}'):
        os.makedirs(f'{path_}/{dssp_folder}', exist_ok=True)
    if not os.path.exists(f'{path_}/{dssp_csv_folder}'):
        os.makedirs(f'{path_}/{dssp_csv_folder}', exist_ok=True)
    fix_pdb(path, f'{path_}/{fix_pdb_folder}')
    dssp_generation(f'{path_}/{fix_pdb_folder}', f'{path_}/{dssp_folder}')
    dsspcsv(f'{path_}/{dssp_folder}', f'{path_}/{dssp_csv_folder}')
    process_protein_files(f'{path_}/{dssp_csv_folder}', name_prefix,
                          f'{path_}/filter_results', start_res, end_res, target_ss, threshold, f'{path_}/fix_pdb')

    return


def check_ss_proportion(csv_path, start_res, end_res, target_ss_list, threshold):
    """
    检查蛋白质指定区域内多个二级结构的综合占比是否超过阈值

    参数:
    csv_path (str): CSV文件路径
    start_res (int): 起始残基序号
    end_res (int): 结束残基序号
    target_ss_list (list): 目标二级结构标识列表 (如['H', 'G'])
    threshold (float): 占比阈值(0-1之间)

    返回:
    tuple: (是否超过阈值, 实际占比, 区域内残基总数)
    """
    # 读取CSV文件
    df = pd.read_csv(csv_path)

    # 筛选指定残基范围内的数据
    region_df = df[(df['Residue_Number'] >= start_res) & (df['Residue_Number'] <= end_res)]

    # 计算区域内残基总数
    total_residues = len(region_df)

    # 如果区域内没有残基，返回False
    if total_residues == 0:
        return False, 0.0, 0

    # 计算目标二级结构的数量（多个标识）
    target_count = region_df[region_df['Secondary_Structure'].isin(target_ss_list)].shape[0]

    # 计算占比
    proportion = target_count / total_residues

    # 判断是否超过阈值
    exceeds_threshold = proportion > threshold

    return exceeds_threshold, proportion, total_residues

def process_protein_files(csv_folder, name_prefix,output_folder, start_res, end_res, target_ss, threshold, fixpdb_folder):
    """
    处理多个蛋白质文件，记录符合条件的蛋白质

    参数:
    file_list (list): 蛋白质CSV文件路径列表
    start_res (int): 起始残基序号
    end_res (int): 结束残基序号
    target_ss (str): 目标二级结构标识
    threshold (float): 占比阈值

    返回:
    list: 符合条件的蛋白质文件路径列表
    """
    positive_hits = []
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    #file_list = os.listdir(csv_folder)
    #folder_path, file_prefix = path.rsplit('/', 1)
    pattern = re.compile(f"^{name_prefix}_\d+\.csv$")
    target_folder = Path(csv_folder)
    matched_files = [
        file for file in target_folder.iterdir()  # 遍历当前文件夹所有条目
        if file.is_file()  # 仅保留文件（过滤文件夹）
           and pattern.match(file.name)  # 正则匹配文件名
    ]
    with open(f'{output_folder}/record.txt', 'w') as f:
        for file_path in matched_files:
            exceeds, prop, count = check_ss_proportion(file_path, start_res, end_res, target_ss, threshold)
            if exceeds:
                print(f'Protein source path: {file_path}')
                print(f'Identification region: {start_res}-{end_res}')
                print(f'Proportion of secondary structure({target_ss}): {prop:.2%}')
                print(f'Exceeds threshold({threshold:.2%}): yes')
                print('')
                f.write(f'Protein source path: {file_path}\n')
                f.write(f'Identification region: {start_res}-{end_res}\n')
                f.write(f'Proportion of secondary structure({target_ss}): {prop:.2%}\n')
                f.write(f'Exceeds threshold({threshold:.2%}): yes\n')
                f.write(f'\n')

                positive_hits.append(file_path)
            else:
                print(f'Protein source path: {file_path}')
                print(f'Identification region: {start_res}-{end_res}')
                print(f'Proportion of secondary structure({target_ss}): {prop:.2%}')
                print(f'Exceeds threshold({threshold:.2%}): no')
                print('')
                f.write(f'Protein source path: {file_path}\n')
                f.write(f'Identification region: {start_res}-{end_res}\n')
                f.write(f'Proportion of secondary structure({target_ss}): {prop:.2%}\n')
                f.write(f'Exceeds threshold({threshold:.2%}): no\n')
                f.write(f'\n')

    for file_path in positive_hits:
        filename = str(file_path).rsplit('/',1)[1]
        file_name = filename.rsplit('.',1)[0]
        shutil.copy(f"{fixpdb_folder}/{file_name}.pdb", f'{output_folder}/{file_name}.pdb')
    return



if __name__ == "__main__":
    args, unknown = parse_args()
    run_rfdiffusion(args, unknown)