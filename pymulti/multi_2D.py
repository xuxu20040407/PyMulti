import subprocess
import os
import shutil
import numpy as np


def init2D(case_dir,program: str):
    pwd = os.getcwd()
    if case_dir==None:
        case_dir = pwd
    else:
        case_dir=os.path.join(pwd, case_dir)

    # 检查case文件夹是否存在，如果不存在则创建
    if not os.path.exists(case_dir):
        os.makedirs(case_dir)

    # 定义新文件夹的路径
    new_dir = os.path.join(case_dir, program)

    # 创建新文件夹
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    # 定义source文件夹的路径
    source_dir = os.path.join(pwd, 'source/2D')

    # 定义要复制的文件列表
    files_to_copy = ['RUN', 'User.r','DEPENDENCES','FILELIST']

    # 检查source文件夹是否存在
    if not os.path.exists(source_dir):
        print(f"错误: 'source'文件夹不存在在路径: {source_dir}")
        return


    # 复制文件到新文件夹
    for file_name in files_to_copy:
        file_path = os.path.join(source_dir, file_name)
        if os.path.exists(file_path):
            shutil.copy(file_path, new_dir)
        else:
            print(f"警告: 文件'{file_name}'在'{source_dir}'中未找到.")

    return source_dir,new_dir


def run_command_2D(new_dir):
    # 构建命令字符串
    command = f"cd {new_dir};chmod 755 ./RUN; ./RUN"

    # 运行命令
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

        # 打印命令的输出
        print(f"命令 '{command}' 的输出:")
        print(result.stdout)

        # 如果有错误输出，也打印出来
        if result.stderr:
            print(f"命令 '{command}' 的错误输出:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        # 如果命令返回非零退出状态，打印错误信息
        print(f"命令 '{command}' 执行失败，返回码：{e.returncode}")
        print(e.stderr)



def data2D_process(program_name,index,stacked_grid):
    pwd=os.getcwd()
    data_dir = os.path.join(pwd,program_name)
    all_data = np.zeros((index, 12))

    for i in range(index):
        folder_path = os.path.join(data_dir, str(i))
        file_path = os.path.join(folder_path, f"fit_{i}.dat")

        if os.path.exists(file_path):
            with open(file_path, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    data = np.fromstring(line, sep=' ')
                    all_data[i, 3:12] = data
    all_data[:, 0] = stacked_grid[:,0]
    all_data[:, 1] = stacked_grid[:,1]
    all_data[:, 2] = stacked_grid[:,2]
    return all_data


