from matplotlib import pyplot as plt

import pymulti as pm
import numpy as np


program_name = 'Multi-1D-MC'
new_dir = pm.init1D(program_name,'fit')
print('init completed')

def external_target_function(sample):
    pm.run_delete_1D(new_dir)
    pm.generate_input_data1D(new_dir, 0,sample)
    pm.run_command_1D(new_dir, 0)
    all_data=pm.data1D_process_fit(program_name,1)
    return all_data[0][2]

laser=[0.000000,
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

min_values = [x * 0.9 for x in laser]
max_values = [x * 1.1 for x in laser]
n_samples = 10

sampler = pm.MC_Sampler(min_values, max_values,optimization_direction='max',T=15)

for i in range(n_samples):
    print("index:"+str(i))
    sampler.MC_sampling(external_target_function)
print('completed')

# 提取样本和对应的函数值
params_samples = np.array(sampler.samples)
values_samples = np.array(sampler.sample_values)
accept_samples = sampler.accepted_samples.copy()
accept_samples.insert(0, 1)
accept_samples = np.array(accept_samples)

np.save('params_samples.npy',params_samples)
np.save('values_samples.npy',values_samples)
np.save('accept_samples.npy',accept_samples)


# 绘制参数的直方图
accepted_params_samples = params_samples[accept_samples == 1]
fig, axs = plt.subplots(2, 1, figsize=(8, 12))
axs[0].hist(accepted_params_samples[:, 0], bins=30, alpha=0.6, color='g')
axs[0].set_title("Histogram of r1")
axs[0].set_xlabel("x")
axs[0].set_ylabel("Density")

axs[1].hist(accepted_params_samples[:, 1], bins=30, alpha=0.6, color='g')
axs[1].set_title("Histogram of r2")
axs[1].set_xlabel("y")
axs[1].set_ylabel("Density")

plt.tight_layout()
plt.savefig('histograms_r2.png')


plt.figure(figsize=(10, 6))
accepted_indices = [i for i, accepted in enumerate(accept_samples) if accepted == 1]
accepted_values = [values_samples[i] for i in accepted_indices]
plt.plot(accepted_indices, accepted_values, 'go', label='Accepted')
rejected_indices = [i for i, accepted in enumerate(accept_samples) if accepted == 0]
rejected_values = [values_samples[i] for i in rejected_indices]
plt.plot(rejected_indices, rejected_values, 'ro', label='Rejected')

plt.title("Optimization Process")
plt.xlabel("Sample Index")
plt.ylabel("Function Value")
plt.legend()
plt.savefig('Optimization2.png')

plt.savefig('Optimization2.png')
plt.show()