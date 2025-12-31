'''
Converts tblout outputs from hmmer to uniprot ids.
'''
import argparse
import subprocess
import sys
import threading
import csv


##input_file是原始蛋白质序列文件，output_id_file是符合阈值的同源序列文件（只有id号），output_file是搜索后的同源序列文件（包含序列），database是蛋白质数据库
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', '-i', type=str)
    parser.add_argument('output_id_file', '-id', type=str)
    parser.add_argument('output_file', '-o', type=str)
    parser.add_argument('database', '-d', type=str)
    return parser.parse_args()


##从tbl文件中提取符合阈值的同源序列 ID，保存到target_ids.txt文件中
def catch_tblid(input_file, output_file):
    with open(input_file) as tbl_file:
        with open(output_file, 'w+') as id_list:
            for line in tbl_file.readlines():
                # skip comments
                if line[0] == '#' or len(line) == 0:
                    continue
                id_list.write(line.split(' ')[0])
                id_list.write('\n')
    return

#seqdb是蛋白质数据库，targetidfile存放需要搜索的蛋白质id文件，output_fasta是输出文件
def catch_sequences(targetidfile, seqdb, output_fasta):
    def run_esl(command):
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
            raise RuntimeError(f"ESL execution failed，exit code: {process.returncode}")
        return
    print('Now create the protein database index ... ')
    command = f'esl-sfetch --index {seqdb}'
    run_esl(command)
    print("[esl-sfetch] fetching sequences from the database ...")
    command = f'esl-sfetch -o {output_fasta} -f {seqdb} {targetidfile}'
    run_esl(command)
    return


def main():
    args = parse_args()
    catch_tblid(args.input_file, args.output_id_file)
    catch_sequences(args.output_id_file, args.database, args.output_fasta)
    return

if __name__ == "__main__":
    main()