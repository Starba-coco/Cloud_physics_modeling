#%%
import pandas as pd
import matplotlib.pyplot as plt

# 데이터 파일 경로
file_path = '/home/dmafytn89/study/Cloud_physics_modeling/Week6/output.txt'

# 데이터 읽기
data = pd.read_csv(file_path, delim_whitespace=True)
data.columns = ['Time', 'w', 'T', 'p', 'z', 'RH', 'Activated_Drops', 'q', 'total_mass']
# data['Activated_Drops'] = data['Activated_Drops'].fillna(0)

# 필요한 열 추출
time = data['Time']
activated_drops = data['total_mass']
# activated_drops = data['Activated_Drops']
# activated_drops = data['RH']
# activated_drops = data['sum_drop']
# activated_drops = data['q']

# 플롯 그리기
plt.figure(figsize=(10, 6))
plt.plot(time, activated_drops, marker='o', linestyle='-', color='b')
plt.xlabel('Time (s)')
# plt.xlim(150, 400)
# plt.ylim(99.5, 100.5)
# plt.ylim(99, 101)
# plt.ylabel('RH (%)')
plt.ylabel('RH (%)')
# plt.ylabel('sum of drops')
plt.ylabel('Total mass')
# plt.ylabel('Activated drops')
plt.title('Total mass over Time')
plt.grid(True)
plt.show()
# %%
