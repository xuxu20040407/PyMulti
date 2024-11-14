from multiprocessing import Pool
import numpy as np
import pymulti as pm  # 确保这个模块已正确安装和导入

n_splits = 2

sampler = pm.Sampler()
Targets_range_min=[0.15, 0.0015, 0.015]
Targets_range_max=[0.17, 0.003, 0.02]
Targets_grid = sampler.uniform_sampling(Targets_range_max, Targets_range_min, n_splits)

# Laser_power_range_min=[0,0,0,10,100]
# Laser_power_range_max=[2,2,5,40,300]
# Laser_time = [0,2,2,2,2,2,2,0]
# Laser_grid = sampler.uniform_sampling([0.17, 0.003, 0.02], [0.15, 0.0015, 0.015], n_splits)
program_name = 'Multi-2D'

def process_task(index, new_dir):
    pm.generate_input_data2D(new_dir, index, None, Targets_grid[index, 0], Targets_grid[index, 1], Targets_grid[index, 2])
    pm.run_command_2D(new_dir, index)

if __name__ == '__main__':
    new_dir = pm.init2D(program_name,n_splits**3)
    with Pool(processes=10) as pool:
        # 使用 functools.partial 来传递 new_dir 参数
        from functools import partial
        pool.map(partial(process_task, new_dir=new_dir), range(n_splits**3))

    inp_data = pm.data2D_process_inp(program_name, n_splits**3)
    fit_data = pm.data2D_process_fit(program_name, n_splits ** 3)
    np.save('2D_database.npy',  np.concatenate((inp_data, fit_data), axis=1))
    pm.run_delete_2D(new_dir)