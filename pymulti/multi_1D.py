import subprocess
import os
import shutil
import numpy as np


def pre_check(program_name):
    # 检查当前项目下是否有tables，如果没有则迁移过来。
    pwd = os.getcwd()
    table_path = os.path.join(program_name, 'tables')
    source_dir = os.path.join(pwd, 'source', 'tables')
    if not os.path.exists(table_path):
        shutil.copytree(source_dir, table_path)


def init1D(program_name):
    # 初始化项目
    pwd = os.getcwd()
    if program_name is None:
        program_name = pwd
    else:
        program_name = os.path.join(pwd, program_name)
    # 创建项目文件夹
    if not os.path.exists(program_name):
        os.makedirs(program_name)
    # 创建项目数据库子文件夹
    new_dir = os.path.join(program_name, 'database')
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    # 迁移源文件至项目文件夹
    source_dir = os.path.join(pwd, 'source/1D')
    files_to_copy = ['fort.12', 'multi']
    # 检查项目文件夹下是否有tables文件夹
    pre_check(program_name)
    for file_name in files_to_copy:
        file_path = os.path.join(source_dir, file_name)
        shutil.copy(file_path, new_dir)

    return new_dir


def generate_input_data1D(case_dir, index, laser, thick1, thick2, thick3):
    if laser is None:
        laser = [0.000000,
                 0.390000,
                 0.570000,
                 0.200000,
                 0.250000,
                 0.300000,
                 0.440000,
                 0.430000,
                 0.330000,
                 0.650000,
                 0.750000,
                 1.350000,
                 0.090000,
                 0.000000,
                 4.730000,
                 0.000000,
                 0.000000,
                 4.000000,
                 4.420000,
                 4.730000,
                 6.860000,
                 18.04000,
                 22.40000,
                 32.00000,
                 32.00000,
                 0.000000]
    if thick1 is None:
        thick1 = 0.15
    if thick2 is None:
        thick2 = 0.003
    if thick3 is None:
        thick3 = 0.02

    if not os.path.exists(case_dir):
        os.makedirs(case_dir)
    # 在指定文件夹下创建算例的输入文件
    output_filename = os.path.join(case_dir, f"inp_{index}.dat")
    with open(output_filename, "w") as fp_out:
        for data in laser:
            fp_out.write(f"{data:.8f}\n")
        fp_out.write(f"{thick1:.4f}\n{thick2:.6f}\n{thick3:.6f}\n")

def run_command_1D(new_dir, index):
    # 执行模拟命令
    command = (f"cd {new_dir};"
               f"rm fit_{index}.dat;"
               f"chmod 755 ./multi;"
               f" ./multi {index}")

    try:
        result = subprocess.run(command, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

        print(f"命令 '{command}' 的输出:")
        print(result.stdout)

        if result.stderr:
            print(f"命令 '{command}' 的错误输出:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print(f"命令 '{command}' 执行失败，返回码：{e.returncode}")
        print(e.stderr)

def data1D_process_inp(program_name, task_num):
    # 先计算inp文件中输入的数目
    pwd = os.getcwd()
    data_dir = os.path.join(pwd, program_name)
    folder_path = os.path.join(data_dir, 'database')
    file_path = os.path.join(folder_path, "inp_0.dat")
    data = np.loadtxt(file_path)
    fit_size=np.size(data,0)
    # 对项目下的所有fit文件进行数据整理
    all_data = np.zeros((task_num, fit_size))
    for i in range(task_num):
        folder_path = os.path.join(data_dir, 'database')
        file_path = os.path.join(folder_path, f"inp_{i}.dat")
        data = np.loadtxt(file_path)
        all_data[i, :] = data
    return all_data

def data1D_process_fit(program_name, task_num):
    # 先计算fit文件中输出的数目
    pwd = os.getcwd()
    data_dir = os.path.join(pwd, program_name)
    folder_path = os.path.join(data_dir, 'database')
    file_path = os.path.join(folder_path, "fit_0.dat")
    data = np.loadtxt(file_path)
    fit_size=np.size(data,0)
    # 对项目下的所有fit文件进行数据整理
    all_data = np.zeros((task_num, fit_size))
    for i in range(task_num):
        folder_path = os.path.join(data_dir, 'database')
        file_path = os.path.join(folder_path, f"fit_{i}.dat")
        data = np.loadtxt(file_path)
        all_data[i, :] = data
    return all_data


def run_delete_1D(new_dir):
    command = f"cd {new_dir}; rm fit_*.dat; rm block_*; rm inp_*.dat"
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

        print(f"命令 '{command}' 的输出:")
        print(result.stdout)

        if result.stderr:
            print(f"命令 '{command}' 的错误输出:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print(f"命令 '{command}' 执行失败，返回码：{e.returncode}")
        print(e.stderr)

