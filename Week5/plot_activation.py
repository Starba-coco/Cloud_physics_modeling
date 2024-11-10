#%%
import pandas as pd
import matplotlib.pyplot as plt

# CSV 파일 읽기
data = pd.read_csv('/home/dmafytn89/study/Cloud_physics_modeling/Week5/output.txt', sep='\+s')

# 필요한 열 추출
time = data.columns[0]
activated_drops = data.columns[6]

# 플롯 그리기
plt.figure(figsize=(10, 6))
plt.plot(time, activated_drops, marker='o', linestyle='-', color='b')
plt.xlabel('Time (s)')
plt.ylabel('Activated Drops')
plt.title('Activated Drops over Time')
plt.grid(True)
plt.show()
