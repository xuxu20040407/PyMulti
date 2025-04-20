import os
import time

def count_dot10_files(directory):
    count = 0
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.10'):
                count += 1
    return count

folder_path = './Multi-1D/database'
initial_dot10_count = count_dot10_files(folder_path)
total_files = 4**8  # 假设总共有4^8个文件
initial_time = time.time()

print(f"There are {initial_dot10_count} files ending with .10 in the directory.")
print(f"Total files expected: {total_files}")
print(f"Initial percentage progress: {initial_dot10_count/total_files*100:.2f}%")

# 设置等待时间为1分钟，即60秒
wait_time = 60

try:
    while time.time() - initial_time < wait_time:
        dot10_count = count_dot10_files(folder_path)
        percentage_progress = dot10_count / total_files * 100
        print(f"Current percentage progress: {percentage_progress:.2f}%")
        time.sleep(10)  # 等待10秒
finally:
    # 1分钟后，计算总运行时间和估算的剩余时间
    elapsed_time = time.time() - initial_time
    estimated_total_time = (total_files - initial_dot10_count) / (dot10_count - initial_dot10_count) * elapsed_time if dot10_count > initial_dot10_count else float('inf')
    estimated_total_time_min=estimated_total_time/60
    print(f"Elapsed time: {elapsed_time:.2f} seconds")
    print(f"Estimated total time to complete: {estimated_total_time:.2f} seconds")
    print(f"Estimated total time to complete: {estimated_total_time_min:.2f} minutes")

