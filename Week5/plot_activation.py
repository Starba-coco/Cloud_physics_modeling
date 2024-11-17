#%%

import pandas as pd
import matplotlib.pyplot as plt

# 데이터 파일 경로
file_path = '/home/dmafytn89/study/Cloud_physics_modeling/Week5/output.txt'

# 데이터 읽기 (공백을 구분자로 지정)
data = pd.read_csv(file_path, delim_whitespace=True)

data.columns = ['Time', 'w', 'T', 'p', 'z', 'RH', 'Activated_Drops', 'q']

# 필요한 열 추출
time = data['Time']
activated_drops = data['Activated_Drops']

# 플롯 그리기
plt.figure(figsize=(10, 6))
plt.plot(time, activated_drops, marker='o', linestyle='-', color='b')
plt.xlabel('Time (s)')
# plt.xlim(1280, 1310)
plt.ylabel('Activated Drops')
plt.title('Activated Drops over Time')
plt.grid(True)
plt.show()

