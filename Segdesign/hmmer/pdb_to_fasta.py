
import argparse
import os
from typing import Dict, List, Optional, Tuple


def pdb_to_fasta(pdb_file, fasta_file):
    """
    从PDB文件中提取蛋白质序列并生成FASTA文件

    参数:
    pdb_file (str): 输入的PDB文件路径
    fasta_file (str): 输出的FASTA文件路径

    # 自动生成同名FASTA文件
    python pdb2fasta.py input.pdb

    # 指定输出FASTA文件名
    python pdb2fasta.py input.pdb -o output.fasta

    # 处理目录中所有PDB文件，输出到fasta_output目录
    python pdb2fasta.py pdb_directory -b

    # 指定输出目录
    python pdb2fasta.py pdb_directory -b -o fasta_results


    """
    # 标准氨基酸的三字母到单字母映射
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # 特殊氨基酸处理
        'SEC': 'U', 'PYL': 'O', 'MSE': 'M', 'SEP': 'S', 'TPO': 'T',
        'PTR': 'Y', 'XLE': 'L', 'GLX': 'Z', 'ASX': 'B'
    }

    # 存储链序列和已处理的残基
    chains = {}
    processed_residues = set()

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                # 只处理ATOM和HETATM记录
                if line.startswith('ATOM'):
                    # 提取链标识符
                    chain_id = line[21]
                    # 提取残基名称
                    res_name = line[17:20].strip()
                    # 提取残基序号和插入码
                    res_seq = line[22:26].strip()
                    icode = line[26].strip()

                    # 创建唯一标识符 (链ID + 残基序号 + 插入码)
                    residue_id = (chain_id, res_seq, icode)

                    # 跳过已处理的残基
                    if residue_id in processed_residues:
                        continue

                    # 检查是否是标准氨基酸
                    if res_name in aa_dict:
                        aa = aa_dict[res_name]
                    else:
                        # 非标准氨基酸用X表示
                        aa = 'X'

                    # 添加到链序列中
                    if chain_id not in chains:
                        chains[chain_id] = []
                    chains[chain_id].append(aa)

                    # 标记此残基已处理
                    processed_residues.add(residue_id)

    except FileNotFoundError:
        print(f"错误: 文件 '{pdb_file}' 未找到")
        return
    output_folder = fasta_file.rsplit('/', 1)[0]
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    # 写入FASTA文件
    with open(fasta_file, 'w') as out:
        for chain_id, seq in chains.items():
            # 生成FASTA头部
            header = f">{pdb_file.split('/')[-1].split('.')[0]}_{chain_id} from {pdb_file}\n"
            out.write(header)

            # 每行80个字符写入序列
            sequence = ''.join(seq)
            for i in range(0, len(sequence), 80):
                out.write(sequence[i:i + 80] + '\n')

    print(f"成功生成FASTA文件: {fasta_file}")


# 示例用法
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PDB文件转FASTA文件工具')
    parser.add_argument('-i', '--input', help='输入PDB文件路径或批量处理的目录')
    parser.add_argument('-o', '--output', help='输出FASTA文件路径（批量处理时为输出目录）')
    parser.add_argument('-b', '--batch', action='store_true', help='批量处理目录中的所有PDB文件')
    #parser.add_argument('-v', '--verbose', action='store_true', help='显示详细处理信息')
    args = parser.parse_args()

    input_file = args.input  # 输入文件
    output_file = args.output  # 输出的FASTA文件
    batch = args.batch
    if batch:
        if not os.path.exists(output_file):  ##新建文件夹
            os.makedirs(output_file, exist_ok=True)
        filenames = os.listdir(input_file)
        for filename in filenames:
            file_name = filename.split('.')[0]
            pdb_to_fasta(f'{input_file}/{filename}', f'{output_file}/{file_name}.fasta')
    else:
        pdb_to_fasta(input_file, output_file)



