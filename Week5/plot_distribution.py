#%%
import numpy as np
import matplotlib.pyplot as plt

# Read data from file
data = np.loadtxt('/home/dmafytn89/study/Cloud_physics_modeling/Week4/distribution.txt', skiprows=1)

# Extract diameter and y_value
diameter_nm = data[:, 0]
y_value = data[:, 1]

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(diameter_nm, y_value, 'b-', linewidth=2)
plt.xlabel('Diameter (nm)')
plt.ylabel(r'$dN/d\log D$ (cm$^{-3}$)')
plt.title('Particle Size Distribution')
plt.xscale('log')
plt.ylim(0, 18)   # Set y-axis limits from 0 to 18
plt.xlim(1, 1000) # Set x-axis limits from 1 to 1000 nm
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.show()
