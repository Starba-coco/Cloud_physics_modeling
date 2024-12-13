#%%
import numpy as np
import matplotlib.pyplot as plt
import glob

path = '/home/dmafytn89/study/Cloud_physics_modeling/Week6/result/'

plt.figure(figsize=(10, 7))

time_steps = np.arange(204, 221, 1)
colors     = plt.cm.viridis(np.linspace(0, 1, len(time_steps)))

for i, t in enumerate(time_steps):
    filename = f'distribution_{t}s.txt'
    try:
        data        = np.loadtxt(path+filename, skiprows=1)
        diameter_nm = data[:, 0]
        y_value     = data[:, 1]

        plt.plot(diameter_nm, y_value, color=colors[i], label=f'{t}s', linewidth=1)
    except FileNotFoundError:
        print(f'File {filename} not found.')

# Customize the plot
plt.xlabel('Diameter (nm)', fontsize=14)
plt.ylabel(r'$dN/d\log D$ (cm$^{-3}$)', fontsize=14)
plt.title('Particle Size Distribution Over Time', fontsize=16)
plt.xscale('log')
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.legend(title='Time', fontsize=10, title_fontsize=12)
plt.tight_layout()
plt.show()

#%%
'''
for ploting drop distribution
'''

import numpy as np
import matplotlib.pyplot as plt

# 데이터를 로드하고 숫자로 변환
path = '/home/dmafytn89/study/Cloud_physics_modeling/Week6/result/'
time_steps = np.arange(0, 3601, 600)
colors = plt.cm.viridis(np.linspace(0, 1, len(time_steps)))

plt.figure(figsize=(10, 7))

for i, t in enumerate(time_steps):
    filename = f'drop_distribution_{t}s.txt'
    try:
        data = np.loadtxt(path + filename, skiprows=1, dtype=str)  # 데이터를 문자열로 읽음
        diameter_nm = np.array([float(d) for d in data[:, 0]])     # 첫 번째 열을 숫자로 변환
        y_value = []

        for val in data[:, 1]:  # 두 번째 열을 숫자로 변환
            try:
                y_value.append(float(val))  # 숫자로 변환 시도
            except ValueError:
                y_value.append(0.0)  # 변환 실패 시 0으로 설정

        y_value = np.array(y_value)

        # 음수를 0으로 변경
        y_value[y_value < 0] = 0

        plt.plot(diameter_nm, y_value, color=colors[i], label=f'{t}s', linewidth=1)
    except FileNotFoundError:
        print(f'File {filename} not found.')

# 그래프 설정
plt.xscale('log')
plt.xlabel('Diameter (nm)')
plt.ylabel('dN/dlogD (cm^-3/nm)')
plt.title('Drop Size Distribution Over Time')
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()


