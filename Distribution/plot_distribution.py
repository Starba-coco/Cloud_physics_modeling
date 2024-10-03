#%%
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('/home/dmafytn89/study/Cloud_physics_modeling/Distribution/result_distribution.txt', delimiter=',', skiprows=1)

r_center = data[:, 0]
n_bin    = data[:, 1] 

plt.figure(figsize=(8, 6))
plt.plot(r_center, n_bin, marker='o', linestyle='-', color='b')
plt.xlabel('Bin Center (m)')
plt.ylabel('Drop Count (#/m³)')
plt.title('Drop Count Distribution')
plt.grid(True)
plt.show()
