import subprocess
import os
import shutil
import numpy as np

def pre_check(case_dir):
    pwd = os.getcwd()
    table_path = os.path.join(case_dir, 'tables')
    source_dir = os.path.join(pwd, 'source', 'tables')
    if not os.path.exists(table_path):
        shutil.copytree(source_dir, table_path)

def init1D(case_dir,program: str):
    pwd = os.getcwd()
    if case_dir==None:
        case_dir = pwd
    else:
        case_dir=os.path.join(pwd, case_dir)

    if not os.path.exists(case_dir):
        os.makedirs(case_dir)

    new_dir = os.path.join(case_dir, program)

    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    source_dir = os.path.join(pwd, 'source/1D')

    files_to_copy = ['fort.12', 'multi']

    if not os.path.exists(source_dir):
        print(f"错误: 'source'文件夹不存在在路径: {source_dir}")
        return
    pre_check(case_dir)

    for file_name in files_to_copy:
        file_path = os.path.join(source_dir, file_name)
        if os.path.exists(file_path):
            shutil.copy(file_path, new_dir)
        else:
            print(f"警告: 文件'{file_name}'在'{source_dir}'中未找到.")

    return source_dir,new_dir


def generate_input_data(case_dir,index,laser,thick1,thick2,thick3):
    if laser==None:
        laser=[
        0.00, 0.13, 2.85, 3.23, 2.31, 1.88, 1.4, 2.4, 2.16, 2.87, 0.92, 2.06, 0.13,
        0.00, 1.8, 1.4, 1.25, 1.12, 3.59, 4.26, 17.89, 29.81, 94.9, 162.11, 300,
        0.00
    ]
    if thick1==None:
        thick1=0.15
    if thick2==None:
        thick2=0.003
    if thick3==None:
        thick3=0.02

    if not os.path.exists(case_dir):
        os.makedirs(case_dir)

    output_filename = os.path.join(case_dir, f"inp_{index}.dat")

    with open(output_filename, "w") as fp_out:
        for data in laser:
            fp_out.write(f"{data:.8f}\n")
        fp_out.write(f"{thick1:.4f}\n{thick2:.6f}\n{thick3:.6f}\n")

def data1D_process(program_name,index,stacked_grid):
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




def run_command_1D(new_dir,index):
    command = f"cd {new_dir};rm fit_{index}.dat;chmod 755 ./multi; ./multi {index}"

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

