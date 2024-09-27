#%%
import pandas as pd
import matplotlib.pyplot as plt

results_path = '/home/dmafytn89/study/Cloud_physics_modeling/results.txt'
results_df   = pd.read_csv(results_path, skipinitialspace=True)

plt.figure(figsize=(10, 6))
plt.scatter(results_df['Time(s)'], results_df['Temperature(K)'], color='blue', s=10)
plt.title('Temperature over time')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.grid(True)
plt.show()
# %%
