import subprocess
import re
import numpy as np

# 测试的 grid_min 范围和步长
start = 0.1
stop = 10
step = 0.1
Average = 100

def run_cuda_program(grid_min):
    command = ["/content/drive/MyDrive/rodinia_3.1/cuda/nn/nn", "list10k.txt", "-g", str(grid_min), "-q", "-r", "10"]
    result = subprocess.run(command, capture_output=True, text=True)
    
    stdout = result.stdout
    time_match = re.search(r'\b(\d+)\s*microseconds\b', stdout)
    if time_match:
        elapsed_time = int(time_match.group(1))
    else:
        elapsed_time = None  
    return elapsed_time

grid_mins = np.arange(start, stop + step, step)
results_str = ""

for grid_min in grid_mins:
    grid_min = round(grid_min, 1)
    times = 0
    for i in range(Average):
        times += run_cuda_program(grid_min)
    average_time = times / Average 
    results_str += f"Average: {Average}, Grid min: {grid_min}, Execution Time: {average_time} ms\n"
    print(f"Average: {Average}, Grid min: {grid_min}, Execution Time: {average_time} ms")
with open('Analysis_Result.txt', 'w') as file:
    file.write(results_str)
