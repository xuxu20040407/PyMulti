import itertools
import os.path
from multiprocessing import Pool
import numpy as np
import pymulti as pm  # 确保这个模块已正确安装和导入
from functools import partial


n_splits = 2

sampler = pm.Sampler()
# Targets=[0.045,0.010,0.000]
# Targets_min=[0.15, 0.0015, 0.015]
# Targets_max=[0.17, 0.003, 0.02]
# Targets_grid = sampler.uniform_sampling(Targets_max, Targets_min, n_splits)

laser_time = [0.000000,
              0.390000,
              0.770000,
              0.250000,
              0.300000,
              0.440000,
              0.430000,
              0.330000,
              0.650000,
              0.750000,
              1.350000,
              0.090000]
laser_power = [0.000000,
               4.730000,
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
laser_time_int=list(itertools.accumulate(laser_time))
laser_time_max=laser_time_int.copy()
laser_time_max[1]=min([laser_time_int[2],laser_time_int[1]+0.1])
laser_time_max[2]=min([laser_time_int[3],laser_time_int[2]+0.1])
laser_time_min=laser_time_int.copy()
laser_time_min[1]=max([laser_time_int[0],laser_time_int[1]-0.1])
laser_time_min[2]=max([laser_time_int[1],laser_time_int[2]-0.1])


laser_power_max=laser_power.copy()
laser_power_max[1]=laser_power[1]*1.5
laser_power_min=laser_power.copy()
laser_power_min[1]=laser_power[1]*0.5

laser_max=laser_time_max+laser_power_max
laser_min=laser_time_min+laser_power_min

Laser_grid = sampler.uniform_sampling(laser_max, laser_min, n_splits)
Laser_grid[:, 1:12] -= Laser_grid[:, :11]
program_name = 'Multi-1D'


def process_task(index, new_dir,Laser_grid):
    pm.generate_input_data1D(new_dir, index,Laser_grid[index])
    pm.run_command_1D(new_dir, index)


if __name__ == '__main__':
    new_dir = pm.init1D(program_name)
    with Pool(processes=10) as pool:
        pool.map(partial(process_task, new_dir=new_dir,Laser_grid=Laser_grid), range(n_splits ** 3))

    inp_data = pm.data1D_process_inp(program_name, n_splits ** 3)
    fit_data = pm.data1D_process_fit(program_name, n_splits ** 3)
    np.save('1D_database.npy', np.concatenate((inp_data, fit_data), axis=1))
    # pm.run_delete_1D(new_dir)
    processor = pm.DataProcessor()
    processor.read_fort10_file(os.path.join(new_dir,'fort_0.10'))
    V = processor.extract('V')
    print(np.min(V))
