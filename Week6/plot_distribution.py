#%%
import numpy as np
import matplotlib.pyplot as plt
import glob

path = '/home/dmafytn89/study/Cloud_physics_modeling/Week6/result/'

plt.figure(figsize=(10, 7))

time_steps = np.arange(205, 221, 1)
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
