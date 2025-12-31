import os
import subprocess
import shutil
import threading
import pandas as pd
import numpy as np
import argparse
import sys
# 获取当前脚本所在目录 (hmmer目录)
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 将根目录添加到 Python 的模块搜索路径中
sys.path.append(root_dir)
from sto2a2m import sto_to_a2m
from tblout2ids import catch_tblid
from tblout2ids import catch_sequences
from pdb_to_fasta import pdb_to_fasta
from dssp.dssp import run_dssp
from dssp.dsspcsv import parse_dssp
#from
from collections import Counter
import math
from scipy.cluster.hierarchy import linkage, fcluster
#from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings("ignore",
                        category=UserWarning,
                        module="evcouplings.couplings.pairs")

#plt.rcParams['font.sans-serif']=['SimHei'] #替换字体
#plt.rcParams['axes.unicode_minus']=False


##--------------外部输入指令--------------------------
def parse_args():
    parser = argparse.ArgumentParser(description='Protein homology analysis',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_pdb', '-p', type=str, help='Path to pdb file')
    parser.add_argument('--select_chain', '-s', type=str, default='A', help='Choose which chain to analyze. If you select chain A, --select_chain A.')
    parser.add_argument('--output_folder', '-o', type=str, help='Output directory')
    parser.add_argument('--bitscore', '-b', type=float, default=0.3, help='Bitscore threshold')
    parser.add_argument('--n_iter', '-n', type=int, default=5, help='Number of iterations')
    parser.add_argument('--database', '-d', type=str, help='Path to database (UniRef100)')
    parser.add_argument('--cpu', '-c', type=int, help='Number of CPU cores used for search')
    parser.add_argument('--minimum_sequence_coverage', type=int, default=50, help='Minimum sequence coverage (percentage)')
    parser.add_argument('--minimum_column_coverage', type=int, default=70, help='Minimum_column_coverage (percentage)')
    parser.add_argument('--identity', '-i', type=float, default=None, help='Identity threshold')
    parser.add_argument('--method', '-m', type=str, default='neff', help='Methods used to analyze non-conserved regions')
    parser.add_argument('--ignore_gaps', '-g', type=bool, default=True, help='Whether to ignore empty placeholders')
    parser.add_argument('--query_sequence', '-q', type=str, default=None, help='Whether to manually input the reference sequence')
    parser.add_argument('--threshold', '-t', type=float, default=6, help='Criteria for identifying non-conserved regions (threshold settings)')
    parser.add_argument("--final_report_folder", type=str, default=None,
                        help="Folder for storing final hmmer_report.csv (default: same as output_folder)")


    return parser.parse_args()
##
#--------------------------从PDB文件中提取指定链的数据并生成新文件--------------------------------------
def extract_chain_from_pdb(input_file, output_file, chain_id='A'):
    """
    从PDB文件中提取指定链的数据并生成新文件

    参数:
    input_file: 输入PDB文件路径
    chain_id: 要提取的链标识符（如'A', 'B', 'C'）
    output_file: 输出PDB文件路径
    """
    # 读取原始文件
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # 筛选数据
    output_lines = []
    chain_id = chain_id.strip().upper()  # 统一转为大写

    for line in lines:
        # 保留所有非原子记录（文件头、注释等）
        if line.startswith(('ATOM', 'HETATM')):
            # 提取链标识符（PDB格式中第22列）
            # 注意：PDB格式中链ID在第22个字符位置（索引21）
            if len(line) >= 22 and line[21] == chain_id:
                output_lines.append(line)


    # 写入新文件
    with open(output_file, 'w') as f:
        f.write('HEADER\n')
        f.write('MODEL        1\n')
        f.writelines(output_lines)
        f.write("ENDMDL".ljust(80) + '\n')
        f.write("END".ljust(80) + '\n')

    print(f"Successfully extracted chain {chain_id}, the new file has been saved as {output_file}")
    print(f"A total of {len(output_lines)} lines of data")
    return




##---------将fasta文件中需要分析的链提取出来，放入一个新的fasta文件中----------------------------
def select_fasta(fasta_file, output_file, select_chain='A'):
    ls = []
    with open(fasta_file, 'r') as f:
        line = f.readline()
        while line:
            ls.append(line)
            line = f.readline()
    target_char = '>'
    # 列表推导式直接生成索引列表
    result_indices = [idx for idx, s in enumerate(ls) if target_char in s]
    index = ord(select_chain) - ord('A')  # A->0, B->1, C->2...
    if 0 <= index < len(result_indices):
        ndx_chain = result_indices[index]
    else:
        raise IndexError(f"The chain index {select_chain} is out of the fasta range (total number of chains in fasta = {len(result_indices)})")

    output_folder = output_file.rsplit('/', 1)[0]
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)

    if index < len(result_indices) - 1:
        with open(output_file, 'a+') as f:
            f.truncate(0)
            for i in range(ndx_chain, result_indices[index + 1]):
                f.write(ls[i])
    else:
        with open(output_file, 'a+') as f:
            f.truncate(0)
            for i in range(ndx_chain, len(ls)):
                f.write(ls[i])
    print(f'The chain that needs to be analyzed has been successfully extracted. The file is stored in:{output_file}')
    return

##----------------------------------------------------------------------------------------------
#################################################################################################

##---------处理输入路径下的fasta文件，将其分割为单序列的fasta文件，支持批量处理--------------------------
def deal_original_fasta(input_fasta, output_doc):

    ###将可能含有多条序列的单个fasta文件进行处理，使每条序列都生成独立的fasta文件，存放在original_single_sequence文件夹中
    def single_fasta(file,s):
        num = 0
        ls = []
        with open(file, 'r') as f:
            line = f.readline()
            while line:
                ls.append(line)
                line = f.readline()
        target_char = '>'
        # 列表推导式直接生成索引列表
        result_indices = [idx for idx, s in enumerate(ls) if target_char in s]
        if output_doc == '':
            if not os.path.exists('original_single_sequence'):  ##新建文件夹
                os.mkdir('original_single_sequence')
            for i in range(len(result_indices) - 1):
                with open(f'original_single_sequence/sequence{i + 1 + s}.fasta', 'a+') as f:
                    f.truncate(0)
                    for j in range(result_indices[i], result_indices[i + 1]):
                        f.write(ls[j])
                print(f'Single sequence file: sequence{i + 1 + s}.fasta has been generated')
                num += 1
            with open(f'original_single_sequence/sequence{len(result_indices) + s}.fasta', 'a+') as f:
                f.truncate(0)
                for j in range(result_indices[-1], len(ls)):
                    f.write(ls[j])
            print(f'Single sequence file: sequence{len(result_indices) + s}.fasta has been generated')
            num += 1
        else:
            os.makedirs(output_doc, exist_ok=True)
            if not os.path.exists(f'{output_doc}/original_single_sequence'):  ##新建文件夹
                os.mkdir(f'{output_doc}/original_single_sequence')
            for i in range(len(result_indices) - 1):
                with open(f'{output_doc}/original_single_sequence/sequence{i + 1 + s}.fasta', 'a+') as f:
                    f.truncate(0)
                    for j in range(result_indices[i], result_indices[i + 1]):
                        f.write(ls[j])
                print(f'Single sequence file: sequence{i + 1 + s}.fasta has been generated')
                num += 1
            with open(f'{output_doc}/original_single_sequence/sequence{len(result_indices) + s}.fasta', 'a+') as f:
                f.truncate(0)
                for j in range(result_indices[-1], len(ls)):
                    f.write(ls[j])
            print(f'Single sequence file: sequence{len(result_indices) + s}.fasta has been generated.')
            num += 1
        return num

    ###如果输入的不是单个fasta文件的路径，而是存放多个fasta文件的文件夹，则对文件夹下的每个fasta文件都进行处理，使每条序列都生成独立的fasta文件，存放在original_single_sequence文件夹中
    def folder():
        filenames = os.listdir(input_fasta)
        n = 0
        for i in filenames:
            n += single_fasta(f'{input_fasta}/{i}',n)
        return

    in_file = input_fasta.split('/')[-1]
    print('---------Preprocessing: now start processing the original input fasta file.------------------------')
    if '.fa' in in_file:
        single_fasta(input_fasta,0)
    else:
        folder()
    print('Preprocessing completed.\n')
    print('')
    return
##----------------------------------------------------------------------------------------------
#################################################################################################

##---------------------------------批量进行同源搜索----------------------------------------
def run_jackhmmer(command):
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

def jackhmmer(input_fasta, output_doc, n_iter, bitscore, database, cpu):

    ###对单序列的fasta文件进行同源搜索
    def jack_single_fasta(fasta_path):

        fasta_name = fasta_path.split('/')[-1].split('.')[0]
        print(f'Searching for homologous sequences in {fasta_name} ... ')
        with open(fasta_path, 'r') as f:
            text = f.read()
        lines = text.split('\n')
        # 2. 提取除第一行外的所有行（索引从1开始）
        other_lines = lines[1:]
        # 3. 将这些行拼接成一行（直接拼接，去除换行符）
        combined_line = ''.join(other_lines)
        seqlen = len(combined_line)
        bitscore_all = seqlen*bitscore

        if output_doc == '':
            if not os.path.exists('tblout'):  ##新建文件夹
                os.mkdir('tblout')
            if not os.path.exists('alignmentfile'):  ##新建文件夹
                os.mkdir('alignmentfile')
            if not os.path.exists('chkhmm'):  ##新建文件夹
                os.mkdir('chkhmm')
            if not os.path.exists('chkali'):  ##新建文件夹
                os.mkdir('chkali')
            if not os.path.exists('logfile'):  ##新建文件夹
                os.mkdir('logfile')
            command = f'jackhmmer -N {n_iter} \
                    --incT {bitscore_all} --incdomT {bitscore_all} -T {bitscore_all} --domT {bitscore_all} \
                    --popen 0.02 --pextend 0.4 --mx BLOSUM62 --tblout tblout/{fasta_name}.tbl \
                    -A alignmentfile/{fasta_name}.sto --noali --notextw --chkhmm chkhmm/{fasta_name} --chkali chkali/{fasta_name} \
                    --cpu {cpu} -o logfile/{fasta_name}_output.log {fasta_path} {database}'
            try:
                run_jackhmmer(command)
                print("Jackhmmer successfully completed！")
            except Exception as e:
                print(f"\nError: {str(e)}")
        else:
            if not os.path.exists(f'{output_doc}/tblout'):  ##新建文件夹
                os.mkdir(f'{output_doc}/tblout')
            if not os.path.exists(f'{output_doc}/alignmentfile'):  ##新建文件夹
                os.mkdir(f'{output_doc}/alignmentfile')
            if not os.path.exists(f'{output_doc}/chkhmm'):  ##新建文件夹
                os.mkdir(f'{output_doc}/chkhmm')
            if not os.path.exists(f'{output_doc}/chkali'):  ##新建文件夹
                os.mkdir(f'{output_doc}/chkali')
            if not os.path.exists(f'{output_doc}/logfile'):  ##新建文件夹
                os.mkdir(f'{output_doc}/logfile')
            command = f'jackhmmer -N {n_iter} \
                    --incT {bitscore_all} --incdomT {bitscore_all} -T {bitscore_all} --domT {bitscore_all} \
                    --popen 0.02 --pextend 0.4 --mx BLOSUM62 --tblout {output_doc}/tblout/{fasta_name}.tbl \
                    -A {output_doc}/alignmentfile/{fasta_name}.sto --noali \
                    --notextw --chkhmm {output_doc}/chkhmm/{fasta_name} --chkali {output_doc}/chkali/{fasta_name} \
                    --cpu {cpu} -o {output_doc}/logfile/{fasta_name}_output.log {fasta_path} {database}'
            try:
                run_jackhmmer(command)
                print("Jackhmmer successfully completed！")
            except Exception as e:
                print(f"\nError: {str(e)}")
        print(f'The homologous sequence search for {fasta_name} is complete.\n')
        return

    #deal_original_fasta(input_fasta, output_doc)

    print('-----------Homology_search: now start using [jackhmmer] for homology search-----------------')
    if output_doc != '':
        os.makedirs(output_doc, exist_ok=True)
    filenames = os.listdir(input_fasta)
    for filename in filenames:
        jack_single_fasta(f'{input_fasta}/{filename}')

    print('Homology search completed.\n')
    return
##-------------------------------------------------------------------------------------------
################################################################################################

##-----------------------同源搜索得到的序列进行筛选，生成相关文件------------------------------
def read_fasta(file_path):
    """
    读取FASTA格式文件，返回序列列表和头部标识列表

    参数:
        file_path (str): FASTA文件路径

    返回:
        tuple: (sequences, headers) 包含两个列表：
               sequences - 序列字符串列表
               headers - 对应的头部标识列表
    """
    sequences = []
    headers = []
    current_seq = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # 保存前一个序列（如果有）
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []

                # 保存新头部标识（去除'>'）
                headers.append(line[1:])
            else:
                # 累积序列行（忽略空行）
                if line:
                    current_seq.append(line)

        # 添加最后一个序列
        if current_seq:
            sequences.append(''.join(current_seq))

    # 验证数据一致性
    if len(sequences) != len(headers):
        raise ValueError(f"FASTA文件格式错误: 序列数量({len(sequences)})与头部标识数量({len(headers)})不匹配")

    return sequences, headers
def data_processing(output_directory, database, identity, minimum_sequence_coverage=50, minimum_column_coverage=70):

    def catch_one_sequences(input_tbl):
        tbl_name = input_tbl.split('/')[-1].split('.')[0]
        catch_tblid(input_file = input_tbl,
                    output_file = f'{output_directory}/data_processing/tblout2ids/homologous_protein_id_file/{tbl_name}_ids.txt')
        catch_sequences(
            targetidfile=f'{output_directory}/data_processing/tblout2ids/homologous_protein_id_file/{tbl_name}_ids.txt',
            seqdb=database,
            output_fasta=f'{output_directory}/data_processing/tblout2ids/homologous_protein_sequence_file/{tbl_name}_homologous_sequences.fasta')
        return

    def sto_to_a2m_one(target_seq_file,input_sto):
        target_seq_file_name = target_seq_file.split('/')[-1].split('.')[0]
        if not os.path.exists(f'{output_directory}/data_processing/sto2a2m/{target_seq_file_name}'):  ##新建文件夹
            os.makedirs(f'{output_directory}/data_processing/sto2a2m/{target_seq_file_name}', exist_ok=True)
        sto_to_a2m(target_seq_file=target_seq_file,
                   sto_alignment_file=input_sto,
                   output_prefix=f'{output_directory}/data_processing/sto2a2m/{target_seq_file_name}/{target_seq_file_name}',
                   minimum_sequence_coverage=minimum_sequence_coverage,
                   minimum_column_coverage=minimum_column_coverage)
        return

    def a2m2a3m(a2m, a3m):
        command = f'reformat.pl {a2m} {a3m}'
        print('a2m2a3m_command:', command)
        result = subprocess.run(command, stdout=subprocess.PIPE, shell=True,
                                stderr=subprocess.PIPE, universal_newlines=True, check=True)
        if os.path.exists(a3m):
            pass
        else:
            a2m_folder = a2m.rsplit('/',1)[0]
            print(f'a2m2a3m_folder: {a2m_folder}')
            #a3m_folder = a3m.rsplit('/')[0]
            a3mname = a3m.rsplit('/',1)[-1]
            print(f'a3m2a3m: {a3mname}')
            a3m_ = os.path.join(a2m_folder, a3mname)
            print(f'a3m2a3m: {a3m_}')
            if os.path.exists(a3m_):
                shutil.move(a3m_, a3m)
            else:
                raise ValueError(f"reformat.pl指令未正常执行, {a3m}未正常生成")

        return

    def load_data(alignment_a3m, output_identities, identity=None):
        align_seqs, align_headers = read_fasta(alignment_a3m)


        df = pd.DataFrame({'id': align_headers, 'sequence': align_seqs})
        df['length'] = df['sequence'].apply(len)

        id_df = pd.read_csv(output_identities)
        df = pd.merge(df, id_df, on='id', how='left')
        if identity is not None:
            df = df[df['identity_to_query'] >= identity]
        return df

    if not os.path.exists(f'{output_directory}/data_processing/tblout2ids/homologous_protein_id_file'):  ##新建文件夹
        os.makedirs(f'{output_directory}/data_processing/tblout2ids/homologous_protein_id_file', exist_ok=True)
    if not os.path.exists(f'{output_directory}/data_processing/tblout2ids/homologous_protein_sequence_file'):  ##新建文件夹
        os.makedirs(f'{output_directory}/data_processing/tblout2ids/homologous_protein_sequence_file', exist_ok=True)
    if not os.path.exists(f'{output_directory}/data_processing/sto2a2m'):  ##新建文件夹
        os.makedirs(f'{output_directory}/data_processing/sto2a2m', exist_ok=True)
    if not os.path.exists(f'{output_directory}/data_processing/a3m'):  ##新建文件夹
        os.makedirs(f'{output_directory}/data_processing/a3m', exist_ok=True)
    if not os.path.exists(f'{output_directory}/data_processing/similarity_filtering'):  ##新建文件夹
        os.makedirs(f'{output_directory}/data_processing/similarity_filtering', exist_ok=True)

    filenames = os.listdir(f'{output_directory}/target_sequence')
    for filename in filenames:
        filename_ = filename.split('.')[0]
        #catch_one_sequences(input_tbl = f'{output_directory}/jackhmmer_out/tblout/{filename_}.tbl')
        sto_to_a2m_one(target_seq_file= f'{output_directory}/target_sequence/{filename}',
                       input_sto = f'{output_directory}/jackhmmer_out/alignmentfile/{filename_}.sto')
        a2m2a3m(a2m=f'{output_directory}/data_processing/sto2a2m/{filename_}/{filename_}.a2m',
                a3m=f'{output_directory}/data_processing/a3m/{filename_}.a3m')
        df = load_data(alignment_a3m=f'{output_directory}/data_processing/a3m/{filename_}.a3m',
                  output_identities=f'{output_directory}/data_processing/sto2a2m/{filename_}/{filename_}_identities.csv', identity=identity)
        df.to_csv(f'{output_directory}/data_processing/similarity_filtering/{filename_}.csv')
        #print('df:',df)

    return

##---------------------------------------------------------------------------------------
######################################################################################################

##-------------------------------------非保守区域分析---------------------------------------------------
# 定义标准氨基酸集合
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AAS = {'B': {'D', 'N'}, 'Z': {'Q', 'E'}, 'J': {'I', 'L'}}

##处理模糊氨基酸残基
def resolve_ambiguous_aa(char, weights=None):
    """
    处理模糊氨基酸残基

    参数:
        char (str): 氨基酸代码
        weights (list): 权重列表（可选）

    返回:
        dict: 可能的标准氨基酸及其权重分布
    """
    if char in AMBIGUOUS_AAS:
        possibilities = AMBIGUOUS_AAS[char]
        if weights:
            # 如果有权重，平均分配到可能的标准氨基酸
            weight_per_aa = weights / len(possibilities)
            return {aa: weight_per_aa for aa in possibilities}
        else:
            # 没有权重信息时，每个可能性计数为1
            return {aa: 1 / len(possibilities) for aa in possibilities}
    return {char: weights if weights else 1}

##计算MSA中每个位置的保守性分数（支持非标准氨基酸）
def calculate_conservation(msa, method='neff', ignore_gaps=False, handle_ambiguous=True):
    """
    计算MSA中每个位置的保守性分数（支持非标准氨基酸）

    参数:
        msa (list of str): 多序列比对结果
        method (str): 计算方法 ('neff' 或 'entropy')
        ignore_gaps (bool): 是否忽略空位符号('-')
        handle_ambiguous (bool): 是否处理模糊氨基酸

    返回:
        np.array: 每个位置的保守性分数
    """
    n_seq = len(msa)
    if n_seq == 0:
        return np.array([])

    n_pos = len(msa[0])

    # 计算序列权重
    weights = calculate_sequence_weights(msa)
    #print('weights', weights)

    # 初始化结果数组
    conservation_scores = np.zeros(n_pos)
    neff_l = np.zeros(n_pos)
    entropy_l = np.zeros(n_pos)

    # 遍历每个位置
    for pos in range(n_pos):
        # 收集该位置的字符及其权重
        char_weights = Counter()
        total_weight = 0.0

        for i, seq in enumerate(msa):
            char = seq[pos]

            # 跳过空位（如果设置）
            if char == '-' and ignore_gaps:
                continue

            w = weights[i]

            # 处理模糊氨基酸
            if handle_ambiguous and char in AMBIGUOUS_AAS:
                resolved = resolve_ambiguous_aa(char, w)
                for aa, aa_weight in resolved.items():
                    char_weights[aa] += aa_weight
                    total_weight += aa_weight
            # 处理未知氨基酸(X)
            elif char == 'X':
                # 排除未知氨基酸
                continue
            # 处理特殊氨基酸(U,O)和终止密码子(*)
            elif char in {'U', 'O', '*'}:
                char_weights[char] += w
                total_weight += w
            # 标准氨基酸
            elif char in STANDARD_AAS:
                char_weights[char] += w
                total_weight += w
            # 其他字符（如小写字母）
            else:
                # 转换为大写处理
                char_upper = char.upper()
                if char_upper in STANDARD_AAS or char_upper in AMBIGUOUS_AAS:
                    char_weights[char_upper] += w
                    total_weight += w
                else:
                    # 未知字符，排除
                    continue

        # 如果没有有效字符，跳过
        if total_weight <= 0:
            conservation_scores[pos] = 0
            continue

        # 计算加权频率
        freq_sum = 0.0
        entropy = 0.0
        #print('char_weights:',char_weights)
        for char, count in char_weights.items():
            freq = count / total_weight
            freq_sum += freq * freq
            #print('freq', freq)
            # 计算熵（仅当freq>0）
            if freq > 0:
                entropy -= freq * math.log2(freq)
                #print('-freq * math.log2(freq):', -freq * math.log2(freq))

        neff = 1.0 / freq_sum if freq_sum > 0 else 0
        neff_l[pos] = neff
        max_entropy = math.log2(20)  # 20种标准氨基酸
        normalized_entropy = entropy / max_entropy
        entropy_l[pos] = normalized_entropy


        #print('entropy', entropy)
        # 计算Neff（有效序列数）
        if method == 'neff':
            # 避免除以零
            #neff = 1.0 / freq_sum if freq_sum > 0 else 0
            conservation_scores[pos] = neff

        # 计算熵
        elif method == 'entropy':
            # 归一化到0-1范围（最大熵是log2(21) - 考虑20个标准氨基酸+特殊字符）
            #max_entropy = math.log2(20)  # 20种标准氨基酸 + 特殊字符
            #normalized_entropy = entropy / max_entropy
            #print('normalized_entropy', normalized_entropy)
            conservation_scores[pos] = normalized_entropy

    #print('conservation_scores', conservation_scores)
    return conservation_scores, neff_l, entropy_l

##计算序列权重
def calculate_sequence_weights(msa):
    """
    计算序列权重（支持非标准氨基酸）

    参数:
        msa (list of str): 多序列比对结果

    返回:
        list: 每个序列的权重
    """
    n_seq = len(msa)
    if n_seq == 0:
        return []
    if n_seq == 1:
        return [1.0]

    n_pos = len(msa[0])

    # 创建包含有效位置的布尔掩码
    valid_positions = []
    for pos in range(n_pos):
        has_valid_char = False
        for seq in msa:
            char = seq[pos].upper()
            if char != '-' and char in STANDARD_AAS:
                has_valid_char = True
                break
        valid_positions.append(has_valid_char)

    # 计算序列间的距离矩阵
    dist_matrix = np.zeros((n_seq, n_seq))
    for i in range(n_seq):
        for j in range(i + 1, n_seq):
            # 计算两个序列的非空位位置的比例差异
            diff_count = 0
            total_pos = 0
            for pos in range(n_pos):
                if not valid_positions[pos]:
                    continue

                c1 = msa[i][pos].upper()
                c2 = msa[j][pos].upper()

                # 跳过非标准氨基酸
                if c1 not in STANDARD_AAS or c2 not in STANDARD_AAS:
                    continue

                # 如果忽略空位，只比较非空位位置
                if c1 != '-' and c2 != '-':
                    total_pos += 1
                    if c1 != c2:
                        diff_count += 1
                # 如果不忽略空位，比较所有位置
                elif c1 != c2:
                    diff_count += 1

            # 处理没有有效位置的情况
            if total_pos == 0:
                distance = 1.0  # 没有可比较位置，设为最大距离
            else:
                distance = diff_count / total_pos

            dist_matrix[i, j] = distance
            dist_matrix[j, i] = distance

    # 转换为压缩距离矩阵
    condensed_dist = []
    for i in range(n_seq):
        for j in range(i + 1, n_seq):
            condensed_dist.append(dist_matrix[i, j])
    condensed_dist = np.array(condensed_dist)

    if condensed_dist.size > 0:
        # 层次聚类（使用平均链接法）
        Z = linkage(condensed_dist, method='average')

        # 动态确定聚类阈值
        mean_distance = np.mean(condensed_dist)
        t = mean_distance if mean_distance > 0 else 0.5

        # 形成聚类
        clusters = fcluster(Z, t=t, criterion='distance')
    else:
        # 没有有效距离，每个序列单独聚类
        clusters = list(range(1, n_seq + 1))

    # 计算每个聚类的权重
    cluster_weights = {}
    for cluster_id in set(clusters):
        cluster_size = sum(clusters == cluster_id)
        cluster_weights[cluster_id] = 1.0 / cluster_size

    # 分配序列权重
    weights = [cluster_weights[cluster_id] for cluster_id in clusters]

    # 归一化权重，使总和为序列数
    weight_sum = sum(weights)
    if weight_sum > 0:
        weights = [w * n_seq / weight_sum for w in weights]

    return weights

##分析MSA中的非保守区域
def analyze_conservation(df, method='neff', threshold=6.0, ignore_gaps=True):
    """
       分析MSA中的非保守区域
       参数:
           df (DataFrame): load_data()返回的数据框
           method (str): 保守性计算方法 ('neff' 或 'entropy')
           threshold (float): 非保守区域的阈值
           ignore_gaps (bool): 是否忽略空位符号

       返回:
           tuple: (conservation_scores, non_conserved_positions)
               conservation_scores: 每个位置的保守性分数列表
               non_conserved_positions: 非保守位置索引列表
       """
    # 获取比对序列
    sequences = df['sequence'].tolist()

    # 验证所有序列长度相同
    seq_lengths = [len(seq) for seq in sequences]
    if len(set(seq_lengths)) > 1:
        raise ValueError("MSA序列长度不一致")

    # 计算保守性分数
    conservation_scores, neff_l , entropy_l = calculate_conservation(
        sequences,
        method=method,
        ignore_gaps=ignore_gaps
    )

    # 识别非保守位置
    non_conserved_positions = []
    for pos, score in enumerate(conservation_scores):
        # 对于Neff: 分数越高表示越不保守
        # 对于熵: 分数越高表示越不保守
        if method == 'neff' and score > threshold:
            non_conserved_positions.append(pos)
        elif method == 'entropy' and score > threshold:
            non_conserved_positions.append(pos)

    return conservation_scores, non_conserved_positions, neff_l , entropy_l

##生成保守性分析报告
def generate_conservation_report(df, conservation_scores, non_conserved_positions, neff_l , entropy_l, ss8, ss3, method='neff', query_sequence=None):
    """
    生成保守性分析报告（支持非标准氨基酸）

    参数:
        df (DataFrame): load_data()返回的数据框
        conservation_scores (list): 每个位置的保守性分数
        non_conserved_positions (list): 非保守位置索引
        query_sequence (str): 查询序列（可选）

    返回:
        DataFrame: 包含详细位置信息的报告
    """
    if query_sequence is None:
        query_sequence = df.iloc[0]['sequence']

    n_pos = len(conservation_scores)
    report_data = []

    for pos in range(n_pos):
        # 获取该位置的氨基酸分布（包括非标准）
        aa_counts = Counter()
        for seq in df['sequence']:
            aa = seq[pos]
            if aa == '-': #空白符号不计入统计
                continue
            aa_counts[aa] += 1

        # 按频率排序
        sorted_aa = sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)
        total = sum(aa_counts.values())

        # 计算最常见的氨基酸占比
        top_aa = sorted_aa[0][0] if sorted_aa else '-'
        top_count = sorted_aa[0][1] if sorted_aa else 0
        top_freq = top_count / total if total > 0 else 0

        # 计算氨基酸类型数量（包括非标准）
        aa_types = len(aa_counts)

        # 计算熵（包括所有字符）
        entropy = 0.0
        if total > 0:
            for count in aa_counts.values():
                freq = count / total
                if freq > 0:
                    entropy -= freq * math.log2(freq)
            max_entropy = math.log2(21)  # 20种标准氨基酸 + 特殊字符
            entropy = entropy / max_entropy

        # 添加特殊标记为非标准残基
        non_standard = ""
        query_aa = query_sequence[pos]
        if query_aa in AMBIGUOUS_AAS:
            non_standard = f"Ambiguous({AMBIGUOUS_AAS[query_aa]})"
        elif query_aa in {'U', 'O', '*', 'X'}:
            non_standard = {
                'U': 'Selenocysteine',
                'O': 'Pyrrolysine',
                '*': 'Stop codon',
                'X': 'Unknown'
            }.get(query_aa, 'Non-standard')

        # 计算保守性水平
        conservation = '.'  # 默认表示保守
        if pos in non_conserved_positions:
            if method == 'neff' and conservation_scores[pos] > 10:
                conservation = 'H'  # High variability
            else:
                conservation = 'M'  # Medium variability

        # 添加到报告
        report_data.append({
            'position': pos + 1,
            'query_aa': query_aa,
            'non_standard': non_standard,
            f'conservation_score({method})': conservation_scores[pos],
            'conservation_level': conservation,
            'top_aa': top_aa,
            'top_freq': top_freq,
            'num_aa_types': aa_types,
            'neff': neff_l[pos],
            'entropy': entropy_l[pos],
            'sequence_context': get_sequence_context(query_sequence, pos),
            f'is_non_conserved(method={method})': pos in non_conserved_positions,
            'ss_8': ss8[pos],
            'ss_3': ss3[pos],
        })

    return pd.DataFrame(report_data)

##获取序列上下文
def get_sequence_context(sequence, center_pos, window=5):
    """
    获取序列上下文

    参数:
        sequence (str): 查询序列
        center_pos (int): 中心位置（0-based）
        window (int): 每侧取多少个氨基酸

    返回:
        str: 序列上下文
    """
    n = len(sequence)
    start = max(0, center_pos - window)
    end = min(n, center_pos + window + 1)

    context = []
    for i in range(start, end):
        if i == center_pos:
            context.append(f"[{sequence[i]}]")  # 标记中心位置
        else:
            context.append(sequence[i])

    return "".join(context)

def secondary_structure(dssp_file):
    with open(dssp_file, 'r') as f:
        content = f.read()
    df = parse_dssp(content)
    df_np = np.array(df)
    ss8 = df_np[:, 3]
    ss3 = df_np[:, 4]
    return ss8, ss3

def summary(report_path, output_folder, file_name, final_report_folder):
    """
    结合二级结构与同源分析，最终给出推荐的非保守区域
    """

    def get_conservation_level(row):
        neff = row["neff"]
        entropy = row["entropy"]

        # 处理缺失值：直接判定为Uncertain
        if pd.isna(neff) or pd.isna(entropy):
            return "Uncertain"

        # 判断High等级
        if neff <= 2.0 and entropy <= 0.2:
            return "High"
        # 判断Medium等级
        elif 2.1 <= neff <= 5.0 and 0.2 < entropy <= 0.4:
            return "Medium"
        # 判断Low等级
        elif 5.0 <= neff <= 15.0 and 0.4 < entropy <= 0.8:
            return "Low"
        # 不符合以上条件则为Uncertain
        else:
            return "Uncertain"

    #path = 'work_test1/conservative_analysis/Dusp4_A_conservation_report.csv'
    df = pd.read_csv(report_path)

    is_non_conserved = [col for col in df.columns if 'is_non_conserved' in col][0]
    is_non_conserved_data = df[is_non_conserved]

    ss8_data = df['ss_8']  # 倒数第一列
    ss3_data = df['ss_3']
    neff = df['neff']
    entropy = df['entropy']
    df["conservation_level"] = df.apply(get_conservation_level, axis=1)
    conservation_level_list = df["conservation_level"].tolist()

    conservation_score = [col for col in df.columns if 'conservation_score' in col][0]
    conservation_score_data = df[conservation_score]

    condition = (is_non_conserved_data == True) & (ss8_data == 'C')
    matching_rows = df[condition].copy()
    matching_indices = matching_rows.index + 1
    # print('matching_indices:',matching_indices)
    indices = matching_indices.tolist()
    # print('indices:', indices)
    df1 = {
        'position': df['position'].tolist(),
        'residue': df['query_aa'].tolist(),
        conservation_score: conservation_score_data.tolist(),
        'ss_8': ss8_data.tolist(),
        'ss_3': ss3_data.tolist(),
        'is_recommended_non_conserved': condition.tolist()
    }
    df2 = {
        'position': df['position'].tolist(),
        'residue': df['query_aa'].tolist(),
        conservation_score: conservation_score_data.tolist(),
        'ss_8': ss8_data.tolist(),
        'ss_3': ss3_data.tolist(),
        'neff': neff.tolist(),
        'entropy': entropy.tolist(),
        'conservation_level': conservation_level_list,
    }
    txt_path = os.path.join(output_folder,f'{file_name}_Recommended_Design_Area.txt')
    with open(txt_path, 'w') as f:
        f.write(
            'Based on the combined analysis of secondary structure and conservation scores, the recommended modification sites are:\n')
        if len(indices) > 0:
            f.write(str(indices).strip('[').strip(']') + '\n')
        else:
            f.write('None\n')

    df1 = pd.DataFrame(df1)
    csv_path = os.path.join(output_folder,f'{file_name}_conservative_comprehensive_report.csv')
    df1.to_csv(csv_path, index=False)
    df2 = pd.DataFrame(df2)
    df2.to_csv(f'{final_report_folder}/hmmer_report.csv')

    return


def conservation_processing(input_pdb, output_path, final_report_folder, method='neff', threshold=6.0, ignore_gaps=True, query_sequence=None):
    print('Begin conservative analysis')
    if not os.path.exists(f'{output_path}/conservative_analysis'):  ##新建文件夹
        os.makedirs(f'{output_path}/conservative_analysis',
                    exist_ok=True)
    filenames = os.listdir(f'{output_path}/data_processing/similarity_filtering')
    for filename in filenames:
        file_name = filename.split('.')[0]

        #pdb_path = f'{output_path}/target_sequence/{file_name}.pdb'
        dssp_path = f'{output_path}/dssp_file/{file_name}.dssp'
        run_dssp(input_pdb, dssp_path)
        ss8, ss3 = secondary_structure(dssp_path)

        df = pd.read_csv(f'{output_path}/data_processing/similarity_filtering/{file_name}.csv')
        conservation_scores, non_conserved_positions, neff_l , entropy_l = analyze_conservation(df, method=method, threshold=threshold,
                                                                            ignore_gaps=ignore_gaps)
        report_df = generate_conservation_report(df, conservation_scores, non_conserved_positions,
                                                 neff_l=neff_l, entropy_l=entropy_l,ss8=ss8, ss3=ss3,
                                                 method=method,query_sequence=query_sequence)
        report_df.to_csv(f'{output_path}/conservative_analysis/{file_name}_conservation_report.csv')
        non_conserved_positions_ndx = report_df[report_df[f'is_non_conserved(method={method})']]['position'].tolist()
        with open(f'{output_path}/conservative_analysis/{file_name}_non_conserved_ndx.txt', 'w') as f:
            f.write('Non-conserved region index:\n')
            f.write(str(non_conserved_positions_ndx).strip('[').strip(']'))
        summary(
            report_path=f'{output_path}/conservative_analysis/{file_name}_conservation_report.csv',
            output_folder=output_path,
            file_name=file_name,
            final_report_folder=final_report_folder,
        )
    print('Conservative analysis completed.')
    return

##--------------------------------------------------------------------------------------------------------
######################################################################################################

##-----------------------Segdesign接口,已弃用------------------------------------
def run_hmmer(param_hmmer):
    input_fasta = param_hmmer['input_fasta']
    output_folder = param_hmmer['output_folder']
    if 'n_iter' in param_hmmer:
        n_iter = int(param_hmmer['n_iter'])
    else:
        n_iter = 5
    if 'bitscore' in param_hmmer:
        bitscore = float(param_hmmer['bitscore'])
    else:
        bitscore = 0.3
    database = param_hmmer['database']
    cpu = int(param_hmmer['cpu'])
    identity = float(param_hmmer['identity'])
    if 'method' in param_hmmer:
        method = param_hmmer['method']
    else:
        method = 'neff'
    if 'threshold' in param_hmmer:
        threshold = float(param_hmmer['threshold'])
    else:
        threshold = 6.0
    if 'ignore_gaps' in param_hmmer:
        ignore_gaps = param_hmmer['ignore_gaps']
    else:
        ignore_gaps = True
    if 'query_sequence' in param_hmmer:
        query_sequence = param_hmmer['query_sequence']
    else:
        query_sequence = None
    deal_original_fasta(input_fasta = input_fasta, output_doc = output_folder)
    jackhmmer(input_fasta=f'{output_folder}/original_single_sequence',
              output_doc=f'{output_folder}/jackhmmer_out',
              n_iter=n_iter,
              bitscore=bitscore,
              database=database,
              cpu=cpu)
    data_processing(output_directory=output_folder,
                    database=database,
                    identity=identity)
    conservation_processing(output_path=output_folder,
                            method=method,
                            threshold=threshold,
                            ignore_gaps=ignore_gaps,
                            query_sequence=query_sequence)
    return



##-------------------------------------主程序-------------------------------------------------------

def main():
    args = parse_args()
    #database = '/home/xieweilong/molecular_dock/hmmer-3.4/install/uniref50_fasta/uniref50.fasta'
    database = os.path.expanduser(args.database)
    input_pdb = os.path.expanduser(args.input_pdb)
    output_folder = os.path.expanduser(args.output_folder)
    final_report_folder = os.path.expanduser(args.final_report_folder)
    fasta_name = args.input_pdb.rsplit('/',1)[-1].strip('.pdb')
    target_chain_pdb_folder = os.path.join(output_folder,'target_chain_pdb')
    if not os.path.exists(target_chain_pdb_folder):
        os.makedirs(target_chain_pdb_folder, exist_ok=True)
    target_chain_pdb_path = os.path.join(target_chain_pdb_folder, f'{fasta_name}_{args.select_chain}.pdb')

    #fasta_path = f'{output_folder}/fasta/{fasta_name}.fasta'
    target_sequence_path = f'{output_folder}/target_sequence/{fasta_name}_{args.select_chain}.fasta'

    extract_chain_from_pdb(input_file=input_pdb, output_file=target_chain_pdb_path, chain_id=args.select_chain)
    pdb_to_fasta(pdb_file=target_chain_pdb_path,
                 fasta_file=target_sequence_path
                 )
    #select_fasta(fasta_file=fasta_path,output_file=target_sequence_path,select_chain=args.select_chain)
    #deal_original_fasta(args.fasta, args.output)
    jackhmmer(input_fasta=target_sequence_path.rsplit('/',1)[0],
              output_doc=f'{output_folder}/jackhmmer_out',
              n_iter=args.n_iter,
              bitscore=args.bitscore,
              database=database,
              cpu=args.cpu
              )
    data_processing(output_directory=output_folder,
                    database=database,
                    identity=args.identity,
                    minimum_sequence_coverage=args.minimum_sequence_coverage,
                    minimum_column_coverage=args.minimum_column_coverage
                    )
    conservation_processing(input_pdb= target_chain_pdb_path,
                            output_path=output_folder,
                            final_report_folder=final_report_folder,
                            method=args.method,
                            threshold=args.threshold,
                            ignore_gaps=args.ignore_gaps,
                            query_sequence=args.query_sequence
                            )




    return

if __name__ == '__main__':
    main()