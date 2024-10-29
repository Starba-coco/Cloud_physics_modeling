#%%
import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the file
data = pd.read_csv('/home/dmafytn89/study/Cloud_physics_modeling/Week4/output.txt', delim_whitespace=True)

# Extract Time and Activated Drops columns
time            = data['Time(s)']
activated_drops = data['Activated']

# Plot Activated Drops over time
plt.figure(figsize=(10, 6))

plt.plot(time, activated_drops, marker='o', linestyle='-', color='b')
plt.xlabel('Time (s)')
plt.ylabel('Activated Drops')
plt.title('Activated Drops over Time')
plt.grid(True)
plt.show()

# %%
