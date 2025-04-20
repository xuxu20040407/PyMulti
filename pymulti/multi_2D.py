import subprocess
import os
import shutil
import numpy as np

import os
import shutil


def init2D(case_dir, task_num,version='fort'):
    # 初始化项目
    pwd = os.getcwd()
    if case_dir is None:
        case_dir = pwd
    else:
        case_dir = os.path.join(pwd, case_dir)
    # 创建项目文件夹
    if not os.path.exists(case_dir):
        os.makedirs(case_dir)
    # 创建项目数据库子文件夹
    new_dir = os.path.join(case_dir, 'database')
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    # 在program下创建0~task_num-1子文件夹
    task_dirs = []
    for i in range(task_num):
        task_dir = os.path.join(new_dir, str(i))
        if not os.path.exists(task_dir):
            os.makedirs(task_dir)
        task_dirs.append(task_dir)
    # 迁移源文件至每个任务子文件夹
    if version == 'fort':
        source_dir = os.path.join(pwd, 'source/2D')
    elif version == 'fit':
        source_dir = os.path.join(pwd, 'source/2D_fit')
    source_dir = os.path.join(pwd, 'source/2D')
    files_to_copy = ['User.r','DEPENDENCES','FILELIST']
    for file_name in files_to_copy:
        file_path = os.path.join(source_dir, file_name)
        if os.path.exists(file_path):
            for task_dir in task_dirs:
                shutil.copy(file_path, task_dir)

    return new_dir


def generate_input_data2D(case_dir, index, laser=None, thick1=0.045, thick2=0.010, thick3=0.000):
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


    if not os.path.exists(case_dir):
        os.makedirs(case_dir)
        os.makedirs(os.path.join(case_dir,str(index)))
    # 在指定文件夹下创建算例的输入文件
    output_filename = os.path.join(case_dir,str(index), f"inp_{index}.dat")
    with open(output_filename, "w") as fp_out:
        for data in laser:
            fp_out.write(f"{data:.8f}\n")
        fp_out.write(f"{thick1:.4f}\n{thick2:.6f}\n{thick3:.6f}\n")


def run_command_2D(new_dir, index):
    # 构建命令字符串
    index_dir=os.path.join(new_dir,str(index))
    command = (f"cd {index_dir};"
               f"rm fit_{index}.dat;"
               f"x.2.library;"
               f"x.2.build multi2d;"
               f"chmod 755 ./multi2d;"
               f"./multi2d {index}")

    # 运行命令
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

        # # 打印命令的输出
        # print(f"命令 '{command}' 的输出:")
        # print(result.stdout)
        #
        # # 如果有错误输出，也打印出来
        # if result.stderr:
        #     print(f"命令 '{command}' 的错误输出:")
        #     print(result.stderr)

    except subprocess.CalledProcessError as e:
        # 如果命令返回非零退出状态，打印错误信息
        print(f"命令 '{command}' 执行失败，返回码：{e.returncode}")
        print(e.stderr)


def data2D_process_inp(program_name, task_num):
    # 先计算inp文件中输入的数目
    pwd = os.getcwd()
    data_dir = os.path.join(pwd, program_name)
    folder_path = os.path.join(data_dir, 'database','0')
    file_path = os.path.join(folder_path, "inp_0.dat")
    data = np.loadtxt(file_path)
    fit_size=np.size(data,0)
    # 对项目下的所有fit文件进行数据整理
    all_data = np.zeros((task_num, fit_size))
    for i in range(task_num):
        folder_path = os.path.join(data_dir, 'database',str(i))
        file_path = os.path.join(folder_path, f"inp_{i}.dat")
        data = np.loadtxt(file_path)
        all_data[i, :] = data
    return all_data

def data2D_process_fit(program_name, task_num):
    # 先计算fit文件中输出的数目
    pwd = os.getcwd()
    data_dir = os.path.join(pwd, program_name)
    folder_path = os.path.join(data_dir, 'database','0')
    file_path = os.path.join(folder_path, "fit_0.dat")
    data = np.loadtxt(file_path)
    fit_size=np.size(data,0)
    # 对项目下的所有fit文件进行数据整理
    all_data = np.zeros((task_num, fit_size))
    for i in range(task_num):
        folder_path = os.path.join(data_dir, 'database',str(i))
        file_path = os.path.join(folder_path, f"fit_{i}.dat")
        data = np.loadtxt(file_path)
        all_data[i, :] = data
    return all_data

