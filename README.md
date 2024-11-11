Pymulti is an interface based on Multi programming, which is designed to facilitate the development of Multi programs in Python.

# File structure
- Pymulti: source code
  -  __init__.py: initialization file
  - Multi _ 1D.py: Methods for 1D
  - Multi _ 2D.py: Methods for 2D
  - Process _ 2D.py: processing method for 2D data
  - Sample _ method. Py: methods for sampling
- Source: The source file for the initialization case
  - 1D: Initialize source files for 1D examples
    - Fort.12: Input Template
    - Multi: Run the program
  - 2D: Initialize the source file for the 2D example
    - Multi2d: Run the program
    - User. R: import templates
  - Tables: files required to run the 1 D program
- Main. Py: sample code

# Function

## init1D(case_dir)
This function creates a new 1D project directory and copies the necessary initialization files from source/1D and source/tables to the database subdirectory of the project.

 **Parameters:**

- Case _ dir (str): Path to the project directory. If not provided, the current working directory will be used.

 **Return**

- New _ dir: Path to database.

For example

```
init1D('Multi-1D')
```

## uniform_sampling(max_values, min_values, n_samples)
This function generates the input data file InP _ { index }.dat for a particular index. The input data includes laser parameters and material thickness.

 **Parameters:**

- Max _ values (list or array): The maximum value for each dimension.
- Min _ values (list or array): The minimum value for each dimension.
- N _ samples (int or list or array): The number of samples per dimension. If it is an integer, all dimensions use the same number of samples; if it is a list or array, it must have the same length as Max _ values and min _ values.

 **Return**
- Sampled _ grid (numpy. Ndarray): 2D array containing sampling grid points.


```
sampler=Sampler()
stacked_grid=sampler.uniform_sampling([0.17,0.003,0.02],[0.15,0.0015,0.015],n_splits)
```
## generate_input_data1D(case_dir, index, laser, thick1, thick2, thick3)
This function generates the input data file InP _ { index }.dat for a particular index. The input data includes laser parameters and material thickness.

 **Parameters:**

- Case _ dir (str): Path to the project directory.
- Index (int): The index of the input file.
- Laser (list, optional): list of laser parameters. The default value is provided.
- Thick1 (float, optional): The thickness of the first layer of material. The default is 0.15.
- Thick2 (float, optional): The thickness of the second layer of material. Default is 0.003.
- Thick3 (float, optional): The thickness of the third layer of material. The default value is 0.02.


```
generate_input_data1D('Multi-1D', 0)
```
The resulting file structure is:
- Multi-1D
  - fort.12
  - multi
  - inp_0.dat

## run_command_1D(new_dir, index)
This function is used to run the calculation task for a specific index. It deletes the old fit _ { index }.dat file, if it exists, and then performs the calculation.

 **Parameters:**

- New _ dir (str): Directory path containing the calculation file.
- Index (int): Compute the index of the task.



```
run_command_1D('Multi-1D', 0)
```
The calculated file structure is:
- Multi-1D
  - fort.12
  - multi
  - inp_0.dat
  - block_0
  - fit_0.dat

## data1D_process(program_name, task_num, stacked_grid)
This function is used to process the fit _ {index}. dat file generated by the calculation and consolidate the data into an array.

 **Parameters:**

- Program _ name (str): The program name used to locate the data directory.
- Task _ num (int): Counts the number of tasks.
- Stacked _ grid (numpy. Ndarray): An array containing grid information.

 **Return**

- All _ data (numpy. Ndarray): An array containing all task data.


```
import numpy as np
stacked_grid = np.array([...])  # 你的网格数据
all_data = data1D_process('Multi-1D', 10, stacked_grid)
```

## run_delete_1D(new_dir)
This function is used to remove all extraneous files after the calculation is complete to keep the directory clean and tidy.

 **Parameters:**

- New _ dir (str): Directory path containing the calculation file.

```
run_delete_1D('Multi-1D')
```
