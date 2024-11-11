import pymulti as pm
from multiprocessing import Pool
import numpy as np

n_splits=2

sampler=pm.Sampler()
stacked_grid=sampler.uniform_sampling([0.17,0.003,0.02],[0.15,0.0015,0.015],n_splits)

program_name='Multi-2D'
def process_task(index):
    
    pm.generate_input_data2D(new_dir, index, None, stacked_grid[index, 0], stacked_grid[index, 1], stacked_grid[index, 2])
    pm.run_command_2D(new_dir, index)

if __name__ == '__main__':
    new_dir = pm.init2D(program_name)
    with Pool(processes=100) as pool:
        pool.map(process_task, range(n_splits**3))
    

    all_data=pm.data2D_process(program_name,n_splits**3,stacked_grid)
    print(all_data)
    np.save('2D_database.npy',all_data)
    pm.run_delete_2D(new_dir)
