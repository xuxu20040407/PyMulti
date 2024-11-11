import pymulti as pm
from multiprocessing import Pool


n_splits=4

sampler=pm.Sampler()
stacked_grid=sampler.uniform_sampling([0.17,0.003,0.02],[0.15,0.0015,0.015],n_splits)

program_name='Multi-1D'
def process_task(index):
    source_dir, new_dir = pm.init1D(program_name, str(index))
    pm.generate_input_data(new_dir, index, None, stacked_grid[index, 0], stacked_grid[index, 1], stacked_grid[index, 2])
    pm.run_command_1D(new_dir, index)

if __name__ == '__main__':
    with Pool(processes=100) as pool:
        pool.map(process_task, range(n_splits**3))

    all_data=pm.data1D_process(program_name,n_splits**3,stacked_grid)