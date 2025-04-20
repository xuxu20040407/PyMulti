from itertools import accumulate
import matplotlib.pyplot as plt
import pymulti as pm
laser_time = [0.000000,
              0.390000,
              0.770000,
              0.250000,
              1.170000,
              0.330000,
              1.400000,
              1.350000,
              0.090000]
laser_power = [0.000000,
               4.730000,
               0.000000,
               4.000000,
               6.860000,
               18.04000,
               32.00000,
               32.00000,
               0.000000]

laser_time_int=list(accumulate(laser_time))
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

sampler = pm.Sampler()
Laser_grid = sampler.uniform_sampling(laser_max, laser_min, 3)
Laser_grid[:, 1:9] -= Laser_grid[:, :8]
# 创建图形
plt.figure(figsize=(10, 5))  # 可以调整图形大小

# 绘制线条
for i in range(21):
    laser_time_int = list(accumulate(Laser_grid[i, 0:9]))
    laser_power = Laser_grid[i, 9:]
    plt.plot(laser_time_int, laser_power, marker='o', linestyle='-', color='b')

# 添加标题和标签
plt.title('Laser Power over Time')
plt.xlabel('Time (s)')
plt.ylabel('Power (W)')

# 显示网格
plt.grid(True)

# 显示图形
plt.show()