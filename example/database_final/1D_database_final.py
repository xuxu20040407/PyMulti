from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import numpy as np
import pymulti as pm
from functools import partial


n_splits = 4

sampler = pm.Sampler()

laser = [0.000000,
	    0.390000,
	    0.770000,
	    0.250000,
	    1.170000,
	    0.330000,
	    1.400000,
	    1.350000,
	    0.090000,
	    0.000000,
	    4.730000,
	    0.000000,
	    4.000000,
	    6.860000,
	    18.04000,
	    32.00000,
	    32.00000,
	    0.000000]

laser_max=[i*1.1 for i in laser]
laser_min=[i*0.9 for i in laser]
laser_max[-3]=32
laser_max[-2]=32
laser_min[-3]=32
laser_min[-2]=32
laser_max[7]=1.350000
laser_max[8]=0.090000
laser_min[7]=1.350000
laser_min[8]=0.090000
laser_max[5]=1.17
laser_min[5]=1.17
laser_max[13]=6.860000
laser_min[13]=6.860000
Laser_grid = sampler.uniform_sampling(laser_max, laser_min, n_splits)

program_name = 'Multi-1D'
print(np.size(Laser_grid, 0))

def process_task(index, new_dir,Laser_grid):
    pm.generate_input_data1D(new_dir, index,Laser_grid[index])
    pm.run_command_1D(new_dir, index)

def thread_task(index, new_dir, Laser_grid):
    with ThreadPoolExecutor(max_workers=4) as executor:
        executor.submit(process_task, index, new_dir, Laser_grid)


if __name__ == '__main__':
    new_dir = pm.init1D(program_name)
    with ProcessPoolExecutor(max_workers=50) as pool:
        pool.map(partial(thread_task, new_dir=new_dir, Laser_grid=Laser_grid), range(n_splits ** 8))


