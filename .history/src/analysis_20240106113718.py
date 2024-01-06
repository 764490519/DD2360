import subprocess
import numpy as np
import matplotlib.pyplot as plt

# 测试的 grid_min 范围和步长
start = 0.1
stop = 100
step = 0.1

def run_cuda_program(grid_min):
    # 替换为您的CUDA程序的路径和必要的命令行参数
    # result = subprocess.run(['./your_cuda_program', str(grid_min)], capture_output=True, text=True)
    result = subprocess.run(['!/content/drive/MyDrive/rodinia_3.1/cuda/nn/nn list10k.txt -q -r 10'], capture_output=True, text=True)
    # 假设程序输出中包含了执行时间
    return float(result.stdout.split()[-1])  # 根据实际输出格式进行调整

grid_mins = np.arange(start, stop + step, step)
execution_times = []

for grid_min in grid_mins:
    time = run_cuda_program(grid_min)
    execution_times.append(time)
    print(f"Grid min: {grid_min}, Execution Time: {time} ms")

# 绘制结果
plt.plot(grid_mins, execution_times)
plt.xlabel('Grid Min Value')
plt.ylabel('Execution Time (ms)')
plt.title('Execution Time vs Grid Min Value')
plt.show()