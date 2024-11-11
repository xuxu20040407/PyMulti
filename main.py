import pymulti as pm

# 调用函数示例
new_dir = pm.init1D('Multi-1D')
pm.generate_input_data1D(new_dir,0, None, None, None, None)
pm.run_command_1D(new_dir,0)

new_dir = pm.init2D('Multi-2D')
pm.generate_input_data2D(new_dir,0, None, None, None, None)
pm.run_command_2D(new_dir,0)
