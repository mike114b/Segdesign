import os
import re
import pandas as pd
import sys
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 将根目录添加到 Python 的模块搜索路径中
sys.path.append(root_dir)
from dssp.dssp import run_dssp
from dssp.dsspcsv import dssp_to_csv


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


def extract_global_score(desc):
    """从FASTA描述信息中提取global_score数值"""
    pattern = re.compile(r'global_score[:=/-]\s*(\d+\.?\d*)', re.IGNORECASE)
    match = pattern.search(desc)
    if match:
        score_str = match.group(1)
        return float(score_str) if '.' in score_str else int(score_str)
    return None


def parse_seq_range(seq_range_str):
    """解析序列区间字符串（如"1-4"），转换为Python切片的start/end索引（0-based）"""
    if not re.match(r'^\d+-\d+$', seq_range_str):
        raise ValueError(f"区间格式错误！请输入如'1-4'/'2-2'的格式，当前输入：{seq_range_str}")

    start_1based, end_1based = map(int, seq_range_str.split('-'))
    if start_1based < 1 or end_1based < 1:
        raise ValueError(f"区间数字必须为正整数！当前输入：{seq_range_str}")
    if start_1based > end_1based:
        raise ValueError(f"起始数字不能大于结束数字！当前输入：{seq_range_str}")

    start_idx = start_1based - 1
    end_idx = end_1based
    return start_idx, end_idx


def parse_uploaded_file(uploaded_file_path):
    """
    解析上传的文件，提取目标名称列表和pLDDT数值
    参数: uploaded_file_path (str): 上传文件的路径（如filter_result.fa.txt）
    返回: tuple: (target_names列表, plddt_dict)
    """
    target_names = []
    plddt_dict = {}
    # 匹配 ">名称, pLDDT=数值" 格式
    pattern = re.compile(r'^>(\S+),\s*pLDDT=(\d+\.?\d*)$')

    with open(uploaded_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):  # 只处理描述行
                match = pattern.match(line)
                if match:
                    name = match.group(1)  # 提取名称（如Dusp4_1_0_0）
                    plddt_val = float(match.group(2))  # 提取pLDDT数值（转换为float）
                    target_names.append(name)
                    plddt_dict[name] = plddt_val
                else:
                    print(f"警告：上传文件中该行格式不匹配，跳过：{line}")

    if not target_names:
        print("警告：上传文件中未提取到有效名称和pLDDT数值")
    return target_names, plddt_dict


def process_fasta_with_filter(folder_path, uploaded_file_path, csv_folder_path, seq_range_str):
    """
    处理fasta文件夹+上传文件筛选，生成五个字典
    参数:
        folder_path (str): fasta文件夹路径
        uploaded_file_path (str): 上传文件路径（如filter_result.fa.txt）
        seq_range_str (str): 序列区间字符串（如"1-4"）
    返回:
        tuple: (筛选后序列字典, 筛选后信息字典, 筛选后global_score字典, 筛选后子序列字典, pLDDT字典)
    """
    # 步骤1：解析上传文件，获取目标名称和pLDDT字典
    target_names, plddt_dict = parse_uploaded_file(uploaded_file_path)
    if not target_names:
        return {}, {}, {}, {}, {}

    # 步骤2：解析序列区间
    try:
        start_idx, end_idx = parse_seq_range(seq_range_str)
    except ValueError as e:
        print(f"区间解析失败：{e}")
        return {}, {}, {}, {}, {}

    # 步骤3：处理fasta文件夹，获取完整字典（同之前逻辑）
    fa_files = [
        f for f in os.listdir(folder_path)
        if f.endswith('.fa') and os.path.isfile(os.path.join(folder_path, f))
    ]
    if not fa_files:
        print("警告：文件夹中未找到.fa文件")
        return {}, {}, {}, {}, {}
    fa_files.sort(key=natural_sort_key)

    def parse_fasta(file_path):
        """解析fasta文件，返回(描述信息, 序列)列表"""
        sequences = []
        current_desc = None
        current_seq = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_desc is not None:
                        sequences.append((current_desc, ''.join(current_seq)))
                    current_desc = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_desc is not None:
                sequences.append((current_desc, ''.join(current_seq)))
        return sequences

    # 生成完整的四个字典
    full_seq_dict = {}
    full_info_dict = {}
    full_global_score_dict = {}
    full_subseq_dict = {}
    s_dict = ss_dict(target_names, csv_folder_path, seq_range_str)
    for file_name in fa_files:
        file_path = os.path.join(folder_path, file_name)
        seq_info_list = parse_fasta(file_path)
        base_name = os.path.splitext(file_name)[0]
        for idx, (desc, seq) in enumerate(seq_info_list):
            new_name = f"{base_name}_{idx}"
            full_seq_dict[new_name] = seq
            full_info_dict[new_name] = desc
            full_global_score_dict[new_name] = extract_global_score(desc)
            full_subseq_dict[new_name] = seq[start_idx:end_idx]

    # 步骤4：按上传文件的目标名称筛选四个字典
    filtered_seq_dict = {name: full_seq_dict[name] for name in target_names if name in full_seq_dict}
    filtered_info_dict = {name: full_info_dict[name] for name in target_names if name in full_info_dict}
    filtered_global_score_dict = {name: full_global_score_dict[name] for name in target_names if
                                  name in full_global_score_dict}
    filtered_subseq_dict = {name: full_subseq_dict[name] for name in target_names if name in full_subseq_dict}

    # 提示未匹配到的名称
    unmatched_names = [name for name in target_names if name not in full_seq_dict]
    if unmatched_names:
        print(f"警告：以下名称在fasta文件夹中未找到匹配项：{unmatched_names}")

    return filtered_seq_dict, filtered_info_dict, filtered_global_score_dict, filtered_subseq_dict, plddt_dict, s_dict


def get_file_dict(folder_path):
    """
    递归读取文件夹下所有文件，生成{无后缀文件名: 绝对路径}的字典
    :param folder_path: 目标文件夹路径
    :return: 文件名-路径字典
    """
    file_dict = {}

    # 验证路径合法性
    if not os.path.exists(folder_path):
        raise FileNotFoundError(f"错误：路径 {folder_path} 不存在")
    if not os.path.isdir(folder_path):
        raise NotADirectoryError(f"错误：{folder_path} 不是有效的文件夹路径")

    # 递归遍历所有文件（包括子文件夹）
    for root, _, files in os.walk(folder_path):
        for file_name in files:
            # 拼接文件绝对路径
            file_abs_path = os.path.abspath(os.path.join(root, file_name))
            # 提取无后缀的文件名（处理多后缀如a.tar.gz时，仅去掉最后一个后缀）
            file_name_no_ext = os.path.splitext(file_name)[0]
            # 存入字典（重复文件名会覆盖）
            file_dict[file_name_no_ext] = file_abs_path

    return file_dict

def ss_dict(target_names, csv_folder_path, seq_range_str):
    s_dict = {}
    start, end = parse_seq_range(seq_range_str)
    print('start:',start)
    print(end)

    file_path_dict = get_file_dict(csv_folder_path)
    for file_name in target_names:
        if file_name in file_path_dict:
            file_path = file_path_dict[file_name]
            df = pd.read_csv(file_path)
            secondary_structure = df["Secondary_Structure"]
            print('secondary_structure:',secondary_structure)

            if start < 0 or end > len(secondary_structure):
                print(f"警告：名称'{file_name}'的范围超出二级结构长度，已跳过")
                continue
            target_struct = ''.join(secondary_structure[start:end].tolist())
            s_dict[file_name] = target_struct
    return s_dict


import csv


def dict_list_to_csv(dict_list, header_list, csv_file_path):
    """
    生成指定样式的CSV文件（支持任意数量属性列/字典）：
    表头：名称       属性0名称    属性1名称    属性2名称    属性3名称 ...
    内容：名称0     属性00       属性10       属性20       属性30
          名称1     属性01       属性11       属性21       属性31
          名称2     属性02       属性12       属性22       属性32

    参数：
    dict_list -- 任意长度的字典列表（每个字典对应1列属性值）
    header_list -- 表头列表（第一个元素为"名称"，后续N个元素对应N个字典的属性名）
    csv_file_path -- 输出CSV文件路径
    """
    # 输入合法性检查（通用适配，不管多少字典都生效）
    if not dict_list:
        raise ValueError("字典列表不能为空！")
    if len(header_list) != len(dict_list) + 1:
        raise ValueError(
            f"表头列表长度应为{len(dict_list) + 1}（1个名称列 + {len(dict_list)}个属性列）！\n"
            f"当前表头长度：{len(header_list)}，字典数量：{len(dict_list)}"
        )
    # 校验所有字典的key是否一致（避免数据错位）
    first_keys = set(dict_list[0].keys())
    for i, d in enumerate(dict_list[1:], 1):
        if set(d.keys()) != first_keys:
            raise ValueError(f"第{i + 1}个字典的key与第一个字典不一致！")

    # 提取行名（名称0、名称1...）并按数字排序
    def sort_key(name):
        num = ''.join([c for c in name if c.isdigit()])
        return int(num) if num else 0

    row_names = sorted(dict_list[0].keys(), key=sort_key)

    # 构建每行数据（自动适配所有字典）
    csv_data = []
    for name in row_names:
        # 一行数据 = [行名] + 该名称在每一个字典中的值（自动遍历所有字典）
        row = [name] + [d[name] for d in dict_list]
        csv_data.append(row)

    # 写入CSV（通用逻辑，无需修改）
    with open(csv_file_path, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(header_list)  # 写入表头
        writer.writerows(csv_data)  # 写入所有行

    print(f"✅ CSV文件已生成：{csv_file_path}")
    print(f"📊 数据维度：{len(row_names)}行（名称0/1/2...） × {len(dict_list)}列（属性列）")

def pdb_to_dssp_csv(pdb_folder, dssp_folder,csv_folder):
    foldernames = os.listdir(pdb_folder)
    for foldername in foldernames:
        folder_path = os.path.join(pdb_folder, foldername)
        dssp_output_folder = os.path.join(dssp_folder, foldername)
        csv_file_path = os.path.join(csv_folder, foldername)
        if not os.path.exists(dssp_output_folder):
            os.makedirs(dssp_output_folder, exist_ok=True)
        if not os.path.exists(csv_file_path):
            os.makedirs(csv_file_path, exist_ok=True)
        for file in os.listdir(folder_path):
            # 拼接完整路径
            file_path = os.path.join(folder_path, file)
            # 判断是否是文件（排除文件夹）
            if os.path.isfile(file_path):
                file_name = os.path.splitext(file)[0]
                dssp_out_file_path = os.path.join(dssp_output_folder, f'{file_name}.dssp')
                csv_out_file_path = os.path.join(csv_file_path, f'{file_name}.csv')
                # 判断文件后缀是否在目标列表中（区分大小写）

                if file_path.endswith('.pdb'):
                    run_dssp(
                        input_path=file_path,
                        output_path=dssp_out_file_path
                    )
                    dssp_to_csv(
                        input_file=dssp_out_file_path,
                        output_file=csv_out_file_path
                    )

    return

def data_organization(
        fasta_folder,
        pdb_folder,
        filter_result_path,
        dssp_folder,
        csv_folder_path,
        seq_range_str
):
    pdb_to_dssp_csv(pdb_folder=pdb_folder, dssp_folder=dssp_folder, csv_folder=csv_folder_path)
    res = process_fasta_with_filter(fasta_folder, filter_result_path, csv_folder_path, seq_range_str)
    filtered_seq, filtered_info, filtered_gs, filtered_subseq, plddt_dict, s_dict = res
    dict_list = [filtered_seq, filtered_subseq, filtered_gs, s_dict, plddt_dict]
    header_list = ['name', 'seq', f'seq({seq_range_str})', 'global_score', f'Secondary_Structure({seq_range_str})',
                   'pLDDT']
    dict_list_to_csv(
        dict_list=dict_list,
        header_list=header_list,
        csv_file_path='out/filter_result.csv')
    return


# 使用示例
if __name__ == "__main__":
    # 配置参数（替换为实际路径）
    fasta_folder = "fasta"  # fasta文件夹路径
    uploaded_file_path = "out/filter_result.fa"  # 上传文件路径
    csv_folder_path = "out/filter_csv"
    seq_range_str = "1-4"  # 序列区间（如"1-4"、"2-2"）

    # 调用函数
    res = process_fasta_with_filter(fasta_folder, uploaded_file_path, csv_folder_path, seq_range_str)
    filtered_seq, filtered_info, filtered_gs, filtered_subseq, plddt_dict, s_dict = res
    print(s_dict)
    # 打印验证结果
    print("=== 上传文件中提取的目标名称 ===")
    print(list(plddt_dict.keys()))

    print(f"\n=== 筛选后-原始序列字典（{len(filtered_seq)}条）===")
    for name, seq in filtered_seq.items():
        print(f"{name}: {seq[:50]}...")

    print(f"\n=== 筛选后-指定区间({seq_range_str})子序列字典 ===")
    for name, sub_seq in filtered_subseq.items():
        print(f"{name}: {sub_seq}")

    print("\n=== 筛选后-信息字典 ===")
    for name, info in filtered_info.items():
        print(f"{name}: {info}")

    print("\n=== 筛选后-global_score字典 ===")
    for name, score in filtered_gs.items():
        print(f"{name}: {score}")

    print("\n=== 第五个-pLDDT字典 ===")
    for name, plddt in plddt_dict.items():
        print(f"{name}: {plddt}")

    print("\n=== 第六个-ss字典 ===")
    for name, ss in s_dict.items():
        print(f"{name}: {ss}")
    #pdb_to_dssp_csv(pdb_folder='out/filter_files', dssp_folder='out/filter_dssp',csv_folder='out/filter_csv')
    dict_list = [filtered_seq, filtered_subseq, filtered_gs, s_dict, plddt_dict]
    header_list = ['name', 'seq', f'seq({seq_range_str})', 'global_score', f'Secondary_Structure({seq_range_str})','pLDDT']
    dict_list_to_csv(
        dict_list=dict_list,
        header_list=header_list,
        csv_file_path='out/filter_result.csv')

    #pdb_to_dssp_csv(pdb_folder='out/filter_files', dssp_folder='out/filter_dssp',csv_folder='out/filter_csv')
    #print(get_file_dict(folder_path='out/filter_csv'))