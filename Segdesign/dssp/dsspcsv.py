import csv
import re

import numpy as np


def parse_dssp(dssp_content):
    """解析DSSP内容，提取残基序号、链标识、氨基酸类型、二级结构代码和描述"""
    data = []
    # 匹配数据行的正则表达式
    pattern = re.compile(r'^\s*(\d+)\s+\d+\s+([A-Z])\s+([A-Z])\s+([A-Z ])')

    # 二级结构代码到描述的映射字典
    ss_dict = {
        'H': 'α-helix',
        'E': 'β-strand',
        'B': 'β-bridge',
        'G': '3/10-helix',
        'I': 'π-helix',
        'P': 'κ-helix',
        'T': 'Turn',
        'S': 'Bend',
        ' ': 'Loop',
        '': 'Loop'  # 处理空字符串的情况
    }
    char_map = {
        'H': 'H',
        'E': 'E',
        'B': 'E',
        'G': 'H',
        'I': 'H',
        'P': 'H',
        'T': 'C',
        'S': 'C',
        'C': 'C',
    }

    for line in dssp_content.split('\n'):
        # 跳过注释行和标题行
        if line.strip().startswith(('#', 'HEADER', 'COMPND', 'SOURCE', '!')):
            continue

        match = pattern.match(line)
        if match:
            res_num = match.group(1).strip()  # 残基序号
            chain = match.group(2).strip()  # 链标识
            amino_acid = match.group(3).strip()  # 氨基酸类型
            ss8_code = match.group(4).strip()  # 二级结构代码
            #print('ss_code:', ss_code)
            # 处理空白的二级结构代码（无规卷曲）
            if ss8_code == '':
                ss8_code = 'C'
                description = 'Loop'
            else:
                # 获取二级结构描述
                description = ss_dict.get(ss8_code, 'Unknown')

            # 如果代码是空格，则使用'C'表示
            if ss8_code == ' ':
                ss8_code = 'C'
            ss3_code = char_map.get(ss8_code, ss8_code)
            data.append([res_num, chain, amino_acid, ss8_code, ss3_code, description])

    return data


def dssp_to_csv(input_file, output_file):
    """读取DSSP文件并输出带有描述的CSV"""
    with open(input_file, 'r') as f:
        content = f.read()

    ss_data = parse_dssp(content)


    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # 添加新列标题
        writer.writerow(['Residue_Number', 'Chain', 'Amino_Acid',
                         'SS_8', 'SS_3', 'SS_Description'])
        writer.writerows(ss_data)

    print(f"Successfully extracted the secondary structure information of {len(ss_data)} residues to {output_file}")


# 使用示例
if __name__ == "__main__":
    input_dssp = "Rad51c.dssp"  # 输入DSSP文件名
    output_csv = "FRad51c_dssp.csv"  # 输出CSV文件名
    dssp_to_csv(input_dssp, output_csv)