#!/usr/bin/env python3
"""
为mpnn_report.csv添加whether_pass列
基于聚类分析输出的代表序列文件确定哪些序列是代表序列
"""

import os
import pandas as pd
from pathlib import Path
import re
import glob
import sys

def extract_representative_indices(result_folder):
    """
    从聚类分析结果文件夹中提取代表序列的index列表
    
    参数:
        result_folder: 聚类结果文件夹路径
        
    返回:
        representative_indices: 代表序列index的集合
    """
    representative_indices = set()
    
    if not os.path.exists(result_folder):
        print(f"警告: 结果文件夹不存在: {result_folder}")
        return representative_indices
    
    # 查找所有FASTA文件（代表序列文件）
    fasta_files = glob.glob(os.path.join(result_folder, "*.fa")) + \
                  glob.glob(os.path.join(result_folder, "*.fasta"))
    
    print(f"📁 在结果文件夹中找到 {len(fasta_files)} 个FASTA文件")
    
    for fasta_file in fasta_files:
        print(f"📄 处理代表序列文件: {os.path.basename(fasta_file)}")
        
        with open(fasta_file, 'r') as f:
            content = f.read()
            lines = content.strip().split('\n')
            
            for line in lines:
                if line.startswith('>'):
                    # 提取序列ID（去掉'>'前缀）
                    seq_id = line[1:].strip()
                    representative_indices.add(seq_id)
                    print(f"  ✅ 添加代表序列: {seq_id}")
    
    print(f"📊 总共找到 {len(representative_indices)} 个代表序列")
    return representative_indices


def add_whether_pass_column(mpnn_report_path, result_folder, output_path=None):
    """
    为mpnn_report.csv添加whether_pass列
    
    参数:
        mpnn_report_path: mpnn_report.csv文件路径
        result_folder: 聚类结果文件夹路径
        output_path: 输出文件路径（默认为覆盖原文件）
        
    返回:
        output_path: 输出文件路径
    """
    print(f"🔄 开始处理: {mpnn_report_path}")
    
    # 读取mpnn_report.csv
    if not os.path.exists(mpnn_report_path):
        raise FileNotFoundError(f"mpnn_report.csv文件不存在: {mpnn_report_path}")
    
    df = pd.read_csv(mpnn_report_path)
    print(f"📊 读取了 {len(df)} 条记录")
    print(f"📋 现有列: {list(df.columns)}")
    
    # 提取代表序列的index列表
    print(f"\n🔍 提取代表序列...")
    representative_indices = extract_representative_indices(result_folder)
    
    # 添加whether_pass列
    print(f"\n➕ 添加whether_pass列...")
    df['whether_pass'] = df['index'].isin(representative_indices)
    
    # 统计结果
    passed_count = df['whether_pass'].sum()
    total_count = len(df)
    
    print(f"📊 统计结果:")
    print(f"  - 总序列数: {total_count}")
    print(f"  - 代表序列数: {passed_count}")
    print(f"  - 非代表序列数: {total_count - passed_count}")
    print(f"  - 代表序列比例: {passed_count/total_count*100:.1f}%")
    
    # 确定输出路径
    if output_path is None:
        output_path = mpnn_report_path
    
    # 保存修改后的文件
    df.to_csv(output_path, index=False)
    print(f"💾 已保存到: {output_path}")
    
    # 显示修改后的前几行
    print(f"\n📄 修改后的前5行数据:")
    print(df[['index', 'whether_pass', 'score']].head())
    
    return output_path


def process_directory_mpnn_reports(mpnn_out_folder):
    """
    批量处理mpnn_out目录下的所有mpnn_report.csv文件
    
    参数:
        mpnn_out_folder: mpnn_out目录路径
    """
    print(f"🔄 批量处理目录: {mpnn_out_folder}")
    
    # 查找所有mpnn_report.csv文件
    report_files = glob.glob(os.path.join(mpnn_out_folder, "**/mpnn_report.csv"), recursive=True)
    
    if not report_files:
        print(f"❌ 在 {mpnn_out_folder} 下未找到mpnn_report.csv文件")
        return
    
    print(f"📁 找到 {len(report_files)} 个mpnn_report.csv文件:")
    for report_file in report_files:
        print(f"  - {report_file}")
    
    # 处理每个文件
    for report_file in report_files:
        print(f"\n" + "="*60)
        
        # 确定对应的result文件夹
        report_dir = os.path.dirname(report_file)
        result_folder = os.path.join(report_dir, "result")
        
        # 检查result文件夹是否存在
        if not os.path.exists(result_folder):
            print(f"⚠️  对应的result文件夹不存在: {result_folder}")
            print(f"   跳过文件: {report_file}")
            continue
        
        try:
            add_whether_pass_column(report_file, result_folder)
            print(f"✅ 成功处理: {report_file}")
        except Exception as e:
            print(f"❌ 处理失败: {report_file}")
            print(f"   错误: {e}")


def test_modification():
    """测试修改功能"""
    print("🧪 测试修改功能...")
    
    # 使用之前创建的测试文件
    # 修正路径：从当前工作目录向上查找
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent  # 从mpnn目录向上两级到项目根目录
    test_dir = project_root / 'test_mpnn_modification'
    mpnn_out_dir = test_dir / 'mpnn_out'
    report_path = mpnn_out_dir / 'mpnn_report.csv'
    result_dir = mpnn_out_dir / 'result'
    
    print(f"🔍 查找测试文件路径:")
    print(f"  - 项目根目录: {project_root}")
    print(f"  - 测试目录: {test_dir}")
    print(f"  - 报告文件: {report_path}")
    print(f"  - 结果目录: {result_dir}")
    
    if not report_path.exists():
        print("❌ 测试文件不存在，请先运行analyze_mpnn_structure.py创建测试数据")
        return False
    
    try:
        output_path = add_whether_pass_column(report_path, result_dir)
        
        # 验证结果
        df_modified = pd.read_csv(output_path)
        if 'whether_pass' in df_modified.columns:
            print("✅ whether_pass列添加成功!")
            
            # 显示统计
            passed = df_modified['whether_pass'].sum()
            total = len(df_modified)
            print(f"📊 测试结果: {passed}/{total} 序列被标记为代表序列")
            
            return True
        else:
            print("❌ whether_pass列添加失败!")
            return False
            
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("=" * 60)
    print("为mpnn_report.csv添加whether_pass列")
    print("=" * 60)
    
    # 测试模式
    if len(sys.argv) == 1:
        print("🧪 运行测试模式...")
        success = test_modification()
        if success:
            print("\n🎉 测试通过!")
        else:
            print("\n❌ 测试失败!")
            sys.exit(1)
    
    # 批量处理模式
    elif len(sys.argv) == 2:
        mpnn_out_folder = sys.argv[1]
        print(f"📁 批量处理模式: {mpnn_out_folder}")
        process_directory_mpnn_reports(mpnn_out_folder)
    
    # 单文件处理模式
    elif len(sys.argv) == 3:
        mpnn_report_path = sys.argv[1]
        result_folder = sys.argv[2]
        print(f"📄 单文件处理模式")
        print(f"  - mpnn_report.csv: {mpnn_report_path}")
        print(f"  - result文件夹: {result_folder}")
        add_whether_pass_column(mpnn_report_path, result_folder)
    
    else:
        print("使用方法:")
        print("  python add_whether_pass_column.py                    # 测试模式")
        print("  python add_whether_pass_column.py <mpnn_out_folder>  # 批量处理模式")
        print("  python add_whether_pass_column.py <report> <result>  # 单文件处理模式")