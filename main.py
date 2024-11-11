import pymulti as pm

# 调用函数示例
source_dir, new_dir = pm.init1D('Multi-1D', '0')
pm.generate_input_data(new_dir, 0, None, None, None, None)
pm.run_command_1D(new_dir, 0)

source_dir, new_dir = pm.init2D('Multi-2D', '0')
pm.run_command_2D(new_dir)
pm.pre