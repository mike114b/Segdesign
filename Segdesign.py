import shutil
import subprocess
import os
import logging
from typing import Dict, Optional, List
import shlex
import argparse
from pathlib import Path
import yaml
import sys
import threading

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(stream=sys.stdout),
        logging.FileHandler('module_runner.log', encoding='utf-8')
    ]
)
logger = logging.getLogger(__name__)

# 配置项（可根据实际情况修改）
CONFIG = {
    "MODULES":{
        'hmmer': {"path":'./Segdesign/hmmer/hmmer.py'},
        'rfdiffusion': {"path":'./Segdesign/rfdiffusion/rf_diffusion.py'},
        'rfdiffusion_report': {"path":'./Segdesign/rfdiffusion/rf_diffusion_report.py'},
        'mpnn': {"path":'./Segdesign/mpnn/mpnn.py'},
        'mpnn_report': {"path":'./Segdesign/mpnn/mpnn_report.py'},
        'esmfold': {"path":'./Segdesign/esmfold/esmfold.py'},
        'esmfold_report': {"path":'./Segdesign/esmfold/esmfold_report.py'},
        'dssp': {"path":'./dssp/dssp.py'},
        'cluster_analysis':{"path":'./Segdesign/mpnn/cluster_analysis.py'},
    },
    "CONFIG_PATH": {
        "MAIN": "./config/config.yaml",
        "SETTING": "./config/setting.yaml"
    }
}



class ModuleRunnerError(Exception):
    """模块运行器自定义异常"""
    pass


def validate_environment(env_name: str) -> bool:
    """验证Conda环境是否存在"""
    conda_info_cmd = [
        f"{CONFIG['MINICONDA_PATH']}/bin/conda",
        "info",
        "--envs"
    ]

    try:
        result = subprocess.run(
            conda_info_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            timeout=30
        )
        # 检查环境是否在输出中（支持完整名称匹配）
        return any(f"*{env_name}" in line or f"  {env_name} " in line for line in result.stdout.splitlines())
    except subprocess.TimeoutExpired:
        logger.warning(f"验证环境 {env_name} 超时")
        return False
    except subprocess.CalledProcessError as e:
        logger.error(f"验证环境失败: {e.stderr}")
        return False


def validate_module(module_name: str) -> str:
    """验证模块是否存在并返回完整路径"""
    if module_name not in CONFIG['MODULES']:
        raise ModuleRunnerError(f"模块 {module_name} 未在配置中定义，可用模块: {list(CONFIG['MODULES'].keys())}")

    module_path = os.path.abspath(CONFIG['MODULES'][module_name]['path'])
    if not os.path.exists(module_path):
        raise ModuleRunnerError(f"模块文件不存在: {module_path}")

    if not os.access(module_path, os.R_OK):
        raise ModuleRunnerError(f"模块文件无读取权限: {module_path}")

    return module_path


def build_command(module_name: str, module_path: str, anaconda_path, env_name: str, custom_args: List[str]) -> str:
    """构建安全的执行命令"""


    # 合并默认参数和自定义参数（自定义参数优先级更高）
    #default_args = MODULE_CONFIG[module_name]["default_args"]
    #final_args = default_args + custom_args

    # 安全转义所有参数，防止命令注入
    escaped_args = [shlex.quote(arg) for arg in custom_args]
    args_str = " ".join(escaped_args)

    # 构建命令（使用set -e确保任一命令失败即退出）
    if anaconda_path is not None:
        anaconda_path = os.path.expanduser(anaconda_path)
        command = f"""
            #!/bin/bash
            set -euo pipefail
            PS1="${{PS1:-}}"
            # 加载conda环境
            if [ -f "{shlex.quote(anaconda_path)}/etc/profile.d/conda.sh" ]; then
                source "{shlex.quote(anaconda_path)}/etc/profile.d/conda.sh"
            elif [ -f "{shlex.quote(anaconda_path)}/bin/activate" ]; then
                source "{shlex.quote(anaconda_path)}/bin/activate"
            else
                echo "找不到conda激活脚本" >&2
                exit 1
            fi

            # 激活环境并运行模块
            conda activate {shlex.quote(env_name)}
            python {shlex.quote(module_path)} {args_str}
            """
    else:
        command = f"""
            # 激活环境并运行模块
            conda run -n {shlex.quote(env_name)} python {shlex.quote(module_path)} {args_str}
            """

    return command
def run_command(command):
    # 创建子进程，捕获标准输出和错误
    print('*'*10)
    print(f"Now starting to execute the command:\n{command}")
    print('*'*10)
    process = subprocess.Popen(
            command,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
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
        raise RuntimeError(f"Command execution failed，exit code: {process.returncode}")
    return


def run_module(
        module_name: str,
        anaconda_path,
        params,
        retry_count: int = 0
) :
    """
    在指定Conda环境中运行模块（支持重试）

    Args:
        module_name: 模块名称
        args: 模块的命令行参数
        retry_count: 当前重试次数

    Returns:
        退出代码（0表示成功）

    Raises:
        ModuleRunnerError: 模块验证或运行失败时抛出
    """
    # 验证模块
    try:
        module_path = validate_module(module_name)
    except ModuleRunnerError as e:
        logger.error(f"模块验证失败: {e}")
        raise

    # 获取环境名称
    env_name = params['env_name']
    logger.info(f"🚀 启动模块: {module_name} (环境: {env_name}, 路径: {module_path})")

    args = [elem for k, v in params['args'].items() for elem in (f'--{k}', str(v))]
    # 构建命令
    command = build_command(
        module_name=module_name,
        module_path=module_path,
        anaconda_path=anaconda_path,
        env_name=env_name,
        custom_args=list(args)
    )

    run_command(command)
    return



def run_module_old(
        module_name: str,
        anaconda_path,
        params,
        retry_count: int = 0
) -> int:
    """
    在指定Conda环境中运行模块（支持重试）

    Args:
        module_name: 模块名称
        args: 模块的命令行参数
        retry_count: 当前重试次数

    Returns:
        退出代码（0表示成功）

    Raises:
        ModuleRunnerError: 模块验证或运行失败时抛出
    """
    # 验证模块
    try:
        module_path = validate_module(module_name)
    except ModuleRunnerError as e:
        logger.error(f"模块验证失败: {e}")
        raise

    # 获取环境名称
    env_name = params['env_name']
    logger.info(f"🚀 启动模块: {module_name} (环境: {env_name}, 路径: {module_path})")

    args = [elem for k, v in params['args'].items() for elem in (f'--{k}', str(v))]
    # 构建命令
    command = build_command(
        module_name=module_name,
        module_path=module_path,
        anaconda_path=os.path.expanduser(anaconda_path),
        env_name=env_name,
        custom_args=list(args)
    )

    try:
        # 执行命令
        result = subprocess.run(
            command,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=CONFIG["COMMAND_TIMEOUT"]
        )

        # 记录输出
        logger.info(f"=== 模块 {module_name} 输出 ===")
        if result.stdout:
            logger.info(result.stdout)
        if result.stderr:
            logger.error(f"模块 {module_name} 错误输出: {result.stderr}")

        logger.info(f"模块 {module_name} 退出代码: {result.returncode}")

        # 重试逻辑
        #if result.returncode != 0 and retry_count < CONFIG["MAX_RETRIES"]:
            #retry_count += 1
            #logger.warning(f"模块 {module_name} 运行失败，将进行第 {retry_count}/{CONFIG['MAX_RETRIES']} 次重试...")
            #return run_module(module_name, *args, retry_count=retry_count)

        return result.returncode

    except subprocess.TimeoutExpired:
        error_msg = f"模块 {module_name} 运行超时（{CONFIG['COMMAND_TIMEOUT']}秒）"
        logger.error(error_msg)
        raise ModuleRunnerError(error_msg) from None
    except subprocess.CalledProcessError as e:
        error_msg = f"模块 {module_name} 运行失败: {e.stderr}"
        logger.error(error_msg)
        raise ModuleRunnerError(error_msg) from e
    except Exception as e:
        error_msg = f"模块 {module_name} 运行异常: {str(e)}"
        logger.error(error_msg, exc_info=True)
        raise ModuleRunnerError(error_msg) from e


def read_yaml_file(yaml_path: str) -> dict:
    """
    读取YAML文件并返回字典格式数据

    Args:
        yaml_path: YAML文件的路径（相对路径或绝对路径）

    Returns:
        解析后的字典数据

    Raises:
        FileNotFoundError: 文件不存在
        yaml.YAMLError: YAML格式错误
        PermissionError: 无文件读取权限
    """
    # 转换为Path对象，方便路径处理
    file_path = Path(yaml_path)

    # 检查文件是否存在
    if not file_path.exists():
        raise FileNotFoundError(f"错误：文件不存在 → {yaml_path}")

    # 检查是否是文件（不是目录）
    if not file_path.is_file():
        raise IsADirectoryError(f"错误：{yaml_path} 是目录，不是文件")

    # 读取并解析YAML文件
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            # yaml.safe_load() 避免执行恶意代码，更安全
            data = yaml.safe_load(f)
        return data or {}
    except PermissionError:
        raise PermissionError(f"错误：无权限读取文件 → {yaml_path}")
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"错误：YAML格式无效 → {e}")
    except Exception as e:
        raise Exception(f"未知错误：{e}")

def merge_configs(config_path: str, setting_path: str) -> dict:
    """
    合并用户配置和系统配置
    
    Args:
        config_path: 用户配置文件路径
        setting_path: 系统配置文件路径
        
    Returns:
        合并后的配置字典
    """
    # 读取配置文件
    user_config = read_yaml_file(config_path)
    setting_config = read_yaml_file(setting_path)
    
    # 合并配置
    merged = {}
    
    # 转换为模块配置
    #merged["modules"] = convert_to_module_config(user_config, setting_config)
    global_parameters = {}
    modules = {}
    project = user_config.get("project", {})
    profile = user_config.get("profile")
    input_pdb =  project.get("input_pdb",'')
    rfdiffusion = user_config.get("rfdiffusion")
    mpnn = user_config.get("mpnn")
    mmseqs = user_config.get("mmseqs")
    esmfold = user_config.get("esmfold")
    output_dir = project.get("output_dir", "./output")

    hmmer_setting = setting_config.get("hmmer", {})  # 无"hmmer"则返回{}
    hmmer_args = hmmer_setting.get("args", {})  # 无"args"则返回{}
    hmmer_user = profile or {}
    hmmer_args.update(hmmer_user)
    hmmer_env = setting_config["environments"].get("hmmer",setting_config["environments"]["main_env"])

    rfdiffusion_setting = setting_config.get("rfdiffusion", {})
    rfdiffusion_args = rfdiffusion_setting.get("args", {})
    rfdiffusion_user = rfdiffusion or {}
    rfdiffusion_args.update(rfdiffusion_user)

    mpnn_setting = setting_config.get("mpnn", {})
    mpnn_args = mpnn_setting.get("args", {})
    mpnn_user = mpnn or {}
    mpnn_args.update(mpnn_user)

    mmseqs_setting = setting_config.get("mmseqs", {})
    mmseqs_args = mmseqs_setting.get("args", {})
    mmseqs_user = mmseqs or {}
    mmseqs_args.update(mmseqs_user)

    esmfold_setting = setting_config.get("esmfold", {})
    esmfold_args = esmfold_setting.get("args", {})
    esmfold_user = esmfold or {}
    esmfold_args.update(esmfold_user)




    # 全局参数配置 (profile)
    if project.get("anaconda_path") is not None:
        global_parameters['anaconda_path'] = project['anaconda_path']
    global_parameters['work_dir'] = output_dir
    merged['global parameters'] = global_parameters

    chain = project.get("chain", "A")

    # hmmer 配置 (profile)
    if profile is not None:
        hmmer_output_folder = os.path.join(output_dir, hmmer_args.get("output_folder", "hmmer_out"))
        hmmer_bitscore = hmmer_args.get("bitscore", 0.3)
        hmmer_n_iter = hmmer_args.get("n_iter", 5)
        hmmer_database = hmmer_args.get("database", "")
        hmmer_cpu = hmmer_args.get("cpu", 10)
        hmmer_minimum_sequence_coverage = hmmer_args.get("minimum_sequence_coverage", 50)
        hmmer_minimum_column_coverage = hmmer_args.get("minimum_column_coverage", 70)
        modules["hmmer"] = {
            "env_name": hmmer_env,
            "args": {
                "input_pdb": input_pdb,
                "select_chain": chain,
                "output_folder": hmmer_output_folder,
                "bitscore": hmmer_bitscore,
                "n_iter": hmmer_n_iter,
                "database": hmmer_database,
                "cpu": hmmer_cpu,
                "minimum_sequence_coverage": hmmer_minimum_sequence_coverage,
                "minimum_column_coverage": hmmer_minimum_column_coverage,
                "final_report_folder": output_dir,  # 新增：最终报告输出到总工作目录
            }
        }


    if project.get("segment") is not None:
        protein_file = os.path.basename(input_pdb)
        protein_name = os.path.splitext(protein_file)[0]

        # rfdiffusion 配置
        if rfdiffusion is not None:
            run_inference_path = rfdiffusion_args["run_inference_path"]
            rfdiffusion_output_folder = os.path.join(output_dir, rfdiffusion_args.get("output_folder","rfdiffusion_out"))
            output_prefix = os.path.join(rfdiffusion_output_folder, f"sample/{protein_name}_{chain}")
            num_designs = rfdiffusion_args.get("num_designs", 10)
            contigs = f"[{project.get('chain', 'A')}1-{project.get('sequence_length', '')}]"
            inpaint_str = f"[{project.get('chain', 'A')}{project.get('segment', '')}]"
            partial_T = rfdiffusion_args["diffuser.partial_T"]
            rfdiffusion_env = setting_config["environments"]["rfdiffusion"]

            modules["rfdiffusion"] = {
                "env_name": rfdiffusion_env,
                "args": {
                    "run_inference_path": run_inference_path,
                    "inference.input_pdb": input_pdb,
                    "inference.output_prefix": output_prefix,
                    "inference.num_designs": num_designs,
                    "contigmap.contigs": contigs,
                    "contigmap.inpaint_str": inpaint_str,
                    "diffuser.partial_T": partial_T
                }
            }
            if rfdiffusion_args.get("contigmap.inpaint_seq") is not None:
                modules["rfdiffusion"]["args"]["contigmap.inpaint_seq"] = rfdiffusion_args.get("contigmap.inpaint_seq")

            # RFdiffusion_report 配置
            rfdiffusion_report_env = setting_config["environments"].get("rfdiffusion_report", setting_config["environments"]["main_env"])
            threshold = rfdiffusion_args.get("threshold", 0.6)
            modules["rfdiffusion_report"] = {
                "env_name": rfdiffusion_report_env,
                "args": {
                    "rfdiffusion_prefix": output_prefix,
                    "inpaint_str": inpaint_str,
                    "threshold": threshold,
                    "final_report_folder": output_dir,  # 新增：最终报告输出到总工作目录
                }

            }

            # 添加结构约束
            select_helix = rfdiffusion_args.get("helix")
            select_strand = rfdiffusion_args.get("strand")
            if select_helix and select_strand is not True:
                modules["rfdiffusion"]["args"]["contigmap.inpaint_str_helix"] = \
                    f"[{project.get('chain', 'A')}{project.get('segment', '')}]"
                modules["rfdiffusion_report"]["args"]['ss'] = f"helix"
            elif select_strand and select_helix is not True:
                modules["rfdiffusion"]["args"]["contigmap.inpaint_str_strand"] = \
                    f"[{project.get('chain', 'A')}{project.get('segment', '')}]"
                modules["rfdiffusion_report"]["args"]['ss'] = "strand"
            else:
                raise ModuleRunnerError(
                    f"Abnormal setting of secondary structure in the design area of module rfdiffusion")

        # mpnn 配置
        if mpnn is not None:
            mpnn_env = setting_config["environments"]["mpnn"]
            parse_multiple_chains_path = mpnn_args["parse_multiple_chains_path"]
            assign_fixed_chains_path = mpnn_args["assign_fixed_chains_path"]
            make_fixed_positions_dict_path = mpnn_args["make_fixed_positions_dict_path"]
            protein_mpnn_run_path = mpnn_args["protein_mpnn_run_path"]
            if mpnn_args.get("pdb_folder") is not None:
                pdb_foler = mpnn_args.get("pdb_folder")
            else:
                pdb_foler = os.path.join(output_dir, f"rfdiffusion_out/filter_results")
            mpnn_output_folder = os.path.join(output_dir, mpnn_args.get("output_folder","mpnn_out"))
            chain_list = project.get("chain", "A")
            position_list =  f"{project.get('chain', 'A')}{project.get('segment', '')}"
            num_seq_per_target = mpnn_args.get("num_seq_per_target", 20)
            sampling_temp = mpnn_args.get("sampling_temp", 0.3)
            seed = mpnn_args.get("seed", 42)

            modules["mpnn"] = {
                "env_name": mpnn_env,
                "args": {
                    "parse_multiple_chains_path": parse_multiple_chains_path,
                    "assign_fixed_chains_path": assign_fixed_chains_path,
                    "make_fixed_positions_dict_path": make_fixed_positions_dict_path,
                    "protein_mpnn_run_path": protein_mpnn_run_path,
                    "pdb_folder": pdb_foler,
                    "output_folder": mpnn_output_folder,
                    "chain_list": chain_list,
                    "position_list": position_list,
                    "num_seq_per_target": num_seq_per_target,
                    "sampling_temp": sampling_temp,
                    "seed": seed,
                    #"top_percent": int(proteinmpnn.get("threshold", 0.9))
                }
            }

            # mpnn_report 配置
            mpnn_report_env = setting_config["environments"].get("mpnn_report",setting_config["environments"]["main_env"])
            seq_folder = os.path.join(mpnn_output_folder, "seqs")
            mpnn_report_output_folder = mpnn_output_folder
            top_percent = mpnn_args.get("top_percent", 0.5)
            rfdiffusion_report_path = mpnn_args.get("rfdiffusion_report_path")

            modules["mpnn_report"] = {
                "env_name": mpnn_report_env,
                "args": {
                    "seq_folder": seq_folder,
                    "output_folder": mpnn_report_output_folder,
                    "top_percent": top_percent,
                    "generate_report": True,  # 添加生成报告标志
                    "final_report_folder": output_dir,  # 新增：最终报告输出到总工作目录
                    "rfdiffusion_report_path": rfdiffusion_report_path
                }
            }
            # 聚类分析配置
            if mmseqs is not None:
                threads = mmseqs_args.get("threads", 8)
                min_seq_id = mmseqs_args.get("min_seq_id")
                cov_mode = mmseqs_args.get("cov_mode", 0)
                coverage = mmseqs_args.get("c", mmseqs_args.get("coverage", 0.8))
                mmseqs_path = mmseqs_args.get("mmseqs_path")
                sensitivity = mmseqs_args.get("s", mmseqs_args.get("sensitivity", 4.0))

                mpnn_report_args_add = {
                    "position_list": position_list,
                    "threads": threads,
                    "min_seq_id": min_seq_id,
                    "cov_mode": cov_mode,
                    "coverage": coverage,
                    "mmseqs_path": mmseqs_path,
                    "sensitivity": sensitivity,
                }
                modules["mpnn_report"]["args"].update(mpnn_report_args_add)
                '''
                modules["mpnn_report"] = {
                    "env_name": mpnn_report_env,
                    "args": {
                        "seq_folder": seq_folder,
                        "output_folder": mpnn_report_output_folder,
                        "top_percent": top_percent,
                        "position_list": position_list,
                        "threads": threads,
                        "min_seq_id": min_seq_id,
                        "cov_mode": cov_mode,
                        "coverage": coverage,
                        "mmseqs_path": mmseqs_path,
                        "sensitivity": sensitivity,

                    }
                }
                '''

        # esmfold 配置
        if esmfold is not None:
            esmfold_env = setting_config["environments"]["esmfold"]
            if esmfold_args.get("input_folder") is not None:
                esmfold_input_folder = esmfold_args.get("input_folder")
            else:
                esmfold_input_folder = os.path.join(output_dir, f"mpnn_out/results")
            esmfold_output_folder = os.path.join(output_dir, esmfold_args.get("output_folder","esmfold_out"))


            modules["esmfold"] = {
                "env_name": esmfold_env,
                "args": {
                    "input_folder": esmfold_input_folder,
                    "output_folder": esmfold_output_folder,
                }
            }

            # esmfold_report 配置
            esmfold_report_env = setting_config["environments"].get("esmfold_report",setting_config["environments"]["main_env"])
            fasta_folder = esmfold_input_folder
            esmfold_folder = esmfold_output_folder
            plddt_threshold = esmfold_args.get("plddt_threshold", 70)
            if esmfold_args.get("original_protein_chain_path") is not None:
                original_protein_chain_path = esmfold_args.get("original_protein_chain_path")
            else:
                chain_folder = os.path.join(output_dir, f"hmmer_out/target_chain_pdb")
                filenames = f"{protein_name}_{chain}.pdb"
                original_protein_chain_path = os.path.join(chain_folder, filenames)

            if esmfold_args.get("seq_range_str") is not None:
                seq_range_str = esmfold_args.get("seq_range_str")
            else:
                seq_range_str = project.get("segment")

            modules["esmfold_report"] = {
                "env_name": esmfold_report_env,
                "args": {
                    "fasta_folder": fasta_folder,
                    "esmfold_folder": esmfold_folder,
                    "plddt_threshold": plddt_threshold,
                    "original_protein_chain_path": original_protein_chain_path,
                    "seq_range_str": seq_range_str,
                }
            }

    # 聚类分析配置
    """
        if project.get("segment") is not None and mmseqs is not None:
        # 动态计算Top百分比文件夹路径
        top_percent_value = mpnn_args.get("top_percent", 0.5)
        top_percent_str = f"{top_percent_value*100:.1f}%"
        
        # 获取mpnn_output_folder，如果不存在则使用默认值
        mpnn_output_folder = os.path.join(output_dir, mpnn_args.get("output_folder", "mpnn_out"))
        top_sequences_folder = os.path.join(mpnn_output_folder, f"top_{top_percent_str}")
        
        # 解析区域位置
        position_range = project.get("segment", "")
        start_pos = int(position_range.split('-')[0]) if '-' in position_range else 1
        end_pos = int(position_range.split('-')[1]) if '-' in position_range else 100
        
        modules["cluster_analysis"] = {
            "env_name": setting_config["environments"].get("cluster_analysis", setting_config["environments"]["main_env"]),
            "args": {
                "input_folder": top_sequences_folder,
                "output_folder": os.path.join(output_dir, "cluster_analysis_out"),
                "start": start_pos,
                "end": end_pos,
                "min_seq_id": mmseqs_args.get("min_seq_id", 0.8),
                "cov_mode": mmseqs_args.get("cov_mode", 0),
                "coverage": mmseqs_args.get("coverage", 0.8),
                "mmseqs_path": mmseqs_args.get("mmseqs_path", "mmseqs"),
                "threads": mmseqs_args.get("threads", 8)
            }
        }

    """

    merged["modules"] = modules
    return merged

def convert_to_module_config(user_config: dict, setting_config: dict) -> dict:
    """
    将用户友好的功能配置转换为模块所需的配置格式
    
    Args:
        user_config: 用户配置
        setting_config: 系统配置
        
    Returns:
        模块配置字典
    """
    modules = {}
    project = user_config.get("project", {})
    profile = user_config.get("profile", {})
    rfdiffusion = user_config.get("rfdiffusion", {})
    proteinmpnn = user_config.get("proteinmpnn", {})
    mmseqs = user_config.get("mmseqs", {})
    esmfold = user_config.get("esmfold", {})
    
    # 输出目录
    output_dir = project.get("output_dir", "./output")
    
    # HMmer 配置 (profile)
    modules["hmmer"] = {
        "env_name": setting_config["environments"]["hmmer"],
        "args": {
            "input_pdb": project.get("input_pdb", ""),
            "select_chain": project.get("chain", ""),
            "output_folder": os.path.join(output_dir, "hmmer_out"),
            "bitscore": profile.get("bitscore", 0.3),
            "n_iter": profile.get("n_iter", 5),
            "database": profile.get("database", ""),
            "cpu": profile.get("cpu", 10),
            "threshold": profile.get("threshold", 0.6)
        }
    }
    # 合并默认参数
    hmmer_config = setting_config.get("hmmer", {})  # 无"hmmer"则返回{}
    hmmer_args = hmmer_config.get("args", {})  # 无"args"则返回{}
    modules["hmmer"]["args"].update(hmmer_args)


    
    # RF Diffusion 配置
    if project.get("segment") is not None:
        modules["rf_diffusion"] = {
            "env_name": setting_config["environments"]["rf_diffusion"],
            "args": {
                "dssp_analyse": ["yes"],
                "threshold": profile.get("threshold", 0.6),
                "run_inference_path": setting_config["rfdiffusion"]["args"]["run_inference_path"],
                "inference.input_pdb": project.get("input_pdb", ""),
                "inference.output_prefix": os.path.join(output_dir, "rfdiffusion_out/sample"),
                "inference.num_designs": rfdiffusion.get("num_designs", 10),
                "contigmap.contigs": [f"{project.get('chain', 'A')}1-{project.get('segment', '').split('-')[1] if '-' in project.get('segment', '') else '100'}"],
                "contigmap.inpaint_str": [f"{project.get('chain', 'A')}{project.get('segment', '')}"],
                "diffuser.partial_T": 50
            }
        }
        
        # 添加结构约束
        if rfdiffusion.get("helix", True):
            modules["rf_diffusion"]["args"]["contigmap.inpaint_str_helix"] = [f"{project.get('chain', 'A')}{project.get('segment', '')}"]
        if rfdiffusion.get("strand", False):
            modules["rf_diffusion"]["args"]["contigmap.inpaint_str_strand"] = [f"{project.get('chain', 'A')}{project.get('segment', '')}"]
    
    # ProteinMPNN 配置
    if project.get("segment") is not None:
        modules["MPNN"] = {
            "env_name": setting_config["environments"]["MPNN"],
            "args": {
                "cluster_analyse": ["yes"],
                "threads": 8,
                "min_seq_id": mmseqs.get("min_seq_id", 0.8),
                "cov_mode": 0,
                "coverage": 0.8,
                "mmseqs_path": "mmseqs",
                "parse_multiple_chains_path": setting_config["MPNN"]["args"]["parse_multiple_chains_path"],
                "assign_fixed_chains_path": setting_config["MPNN"]["args"]["assign_fixed_chains_path"],
                "make_fixed_positions_dict_path": setting_config["MPNN"]["args"]["make_fixed_positions_dict_path"],
                "protein_mpnn_run_path": setting_config["MPNN"]["args"]["protein_mpnn_run_path"],
                "pdb_path": os.path.join(output_dir, "rfdiffusion_out"),
                "output_folder": os.path.join(output_dir, "mpnn_out"),
                "chain_list": project.get("chain", ""),
                "position_list": f"{project.get('chain', 'A')}{project.get('segment', '')}",
                "num_seq_per_target": proteinmpnn.get("num_seq_per_target", 20),
                "sampling_temp": proteinmpnn.get("sampling_temp", 0.3),
                "seed": proteinmpnn.get("seed", 42),
                "top_percent": int(proteinmpnn.get("threshold", 0.9) * 100)
            }
        }
    
    # ESMFold 配置
    if project.get("segment") is not None:
        modules["esmfold"] = {
            "env_name": setting_config["environments"]["esmfold"],
            "args": {
                "input_folder": os.path.join(output_dir, "mpnn_out/top_90.0%"),
                "output_folder": os.path.join(output_dir, "esmfold_out"),
                "plddt_threshold": esmfold.get("plddt_threshold", 70)
            }
        }
    
    
    return modules

def global_work_dir_handling(yaml_data):
    """处理工作目录"""
    work_dir = os.path.expanduser(yaml_data.get('global parameters', {}).get("work_dir", "./output"))
    if not os.path.exists(work_dir):
        os.makedirs(work_dir, exist_ok=True)
    return work_dir







if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SegDesign: 蛋白质设计工具",
        epilog="示例：python Segdesign.py --config ./config/config.yaml --setting ./config/setting.yaml"
    )

    # 添加参数
    parser.add_argument(
        "--config",
        type=str,
        default=CONFIG["CONFIG_PATH"]["MAIN"],
        help="用户配置文件路径（相对路径或绝对路径）"
    )
    parser.add_argument(
        "--setting",
        type=str,
        default=CONFIG["CONFIG_PATH"]["SETTING"],
        help="系统配置文件路径（相对路径或绝对路径）"
    )
    
    args = parser.parse_args()
    
    try:
        # 合并配置
        merged_config = merge_configs(args.config, args.setting)
        print("✅ 配置文件读取成功！")
        print("📊 解析后的数据：")
        print(yaml.dump(merged_config, allow_unicode=True, sort_keys=False))
        
        # 处理工作目录
        output_dir = global_work_dir_handling(merged_config)
        logger.info(f"工作目录: {output_dir}")

        #将config.yaml复制到工作目录下

        shutil.copy(args.config, f"{output_dir}/config.yaml")
        
        # 获取anaconda路径
        anaconda_path = merged_config["global parameters"].get("anaconda_path")
        
        # 运行模块
        for module_name, params in merged_config["modules"].items():
            if module_name in CONFIG['MODULES']:
                try:
                    logger.info(f"正在运行模块: {module_name}")
                    run_module(
                        module_name=module_name,
                        anaconda_path=anaconda_path,
                        params=params
                    )
                    logger.info(f"✅ 模块 {module_name} 运行成功")
                except ModuleRunnerError as e:
                    logger.critical(f"❌ 模块 {module_name} 运行失败: {e}")
                    exit(1)
                except KeyboardInterrupt:
                    logger.info("程序被用户中断")
                    exit(0)
                except Exception as e:
                    logger.critical(f"❌ 模块 {module_name} 未预期的错误: {str(e)}", exc_info=True)
                    exit(1)
        
        logger.info("🎉 所有模块运行完成！")
        
    except Exception as e:
        print(f"❌ 程序执行失败：{e}")
        logger.error(f"程序执行失败: {e}", exc_info=True)
        exit(1)  # 非0退出码表示程序异常


