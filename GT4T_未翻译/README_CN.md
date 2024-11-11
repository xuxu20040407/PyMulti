Pymulti 是基于Multi程序设计的接口，旨在方便Multi程序在Python下的开发。

# 文件结构
- pymulti： 源代码
  - __init__.py： 初始化文件
  - multi_1D.py： 有关1D的方法
  - multi_2D.py： 有关2D的方法
  - process_2D.py： 有关2D数据的处理方法
  - sample_method.py： 有关采样的方法
- source： 初始化算例的源文件
  - 1D： 初始化1D算例的源文件
    - fort.12：输入模板
    - multi：运行程序
  - 2D： 初始化2D算例的源文件
    - multi2d：运行程序
    - User.r： 输入模板
  - tables： 运行1D程序必备的文件
- main.py： 示例代码

# 函数

## init1D(case_dir)
此函数用于创建一个新的1D项目目录，并将必要的初始化文件从source/1D和source/tables复制到项目的database子目录中。

**参数：**

- case_dir (str): 项目目录的路径。如果未提供，将使用当前工作目录。

**返回：**

- new_dir: database的路径。

比如
```
init1D('Multi-1D')
```

## uniform_sampling(max_values, min_values, n_samples)
此函数用于生成特定索引的输入数据文件inp_{index}.dat。输入数据包括激光参数和材料厚度。

**参数：**

- max_values (list or array): 每个维度的最大值。
- min_values (list or array): 每个维度的最小值。
- n_samples (int or list or array): 每个维度的样本数。如果为整数，则所有维度使用相同的样本数；如果为列表或数组，则必须与max_values和min_values具有相同的长度。

**返回：**
- sampled_grid (numpy.ndarray): 包含采样网格点的2D数组。

```
sampler=Sampler()
stacked_grid=sampler.uniform_sampling([0.17,0.003,0.02],[0.15,0.0015,0.015],n_splits)
```
## generate_input_data1D(case_dir, index, laser, thick1, thick2, thick3)
此函数用于生成特定索引的输入数据文件inp_{index}.dat。输入数据包括激光参数和材料厚度。

**参数：**

- case_dir (str): 项目目录的路径。
- index (int): 输入文件的索引。
- laser (list, optional): 激光参数列表。默认值已提供。
- thick1 (float, optional): 第一层材料的厚度。默认值为0.15。
- thick2 (float, optional): 第二层材料的厚度。默认值为0.003。
- thick3 (float, optional): 第三层材料的厚度。默认值为0.02。

```
generate_input_data1D('Multi-1D', 0)
```
生成的文件结构为：
- Multi-1D
  - fort.12
  - multi
  - inp_0.dat

## run_command_1D(new_dir, index)
此函数用于运行特定索引的计算任务。它将删除旧的fit_{index}.dat文件（如果存在），然后执行计算。

**参数：**

- new_dir (str): 包含计算文件的目录路径。
- index (int): 计算任务的索引。


```
run_command_1D('Multi-1D', 0)
```
计算后的文件结构为：
- Multi-1D
  - fort.12
  - multi
  - inp_0.dat
  - block_0
  - fit_0.dat

## data1D_process(program_name, task_num, stacked_grid)
此函数用于处理计算生成的fit_{index}.dat文件，并将数据整合到一个数组中。

**参数：**

- program_name (str): 程序名称，用于定位数据目录。
- task_num (int): 计算任务的数量。
- stacked_grid (numpy.ndarray): 包含网格信息的数组。

**返回：**

- all_data (numpy.ndarray): 包含所有任务数据的数组。

```
import numpy as np
stacked_grid = np.array([...])  # 你的网格数据
all_data = data1D_process('Multi-1D', 10, stacked_grid)
```

## run_delete_1D(new_dir)
此函数用于在计算完成后删除所有无关文件，以保持目录的整洁。

**参数：**

- new_dir (str): 包含计算文件的目录路径。
```
run_delete_1D('Multi-1D')
```
