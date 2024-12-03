#%%
import numpy as np
import matplotlib.pyplot as plt

# 파일 경로
file_path = "/home/dmafytn89/study/Cloud_physics_modeling/Week6/vt_output.txt"

# 데이터 로드
try:
    data = np.loadtxt(file_path, skiprows=1)  # 첫 번째 행은 헤더이므로 건너뜀
    radii = data[:, 0]    # Radius (m)
    vt = data[:, 1]       # Terminal Velocity (cm/s)
except Exception as e:
    print(f"Error loading file: {e}")
    exit()

# 플롯 생성
plt.figure(figsize=(10, 6))
plt.plot(radii, vt, marker='o', linestyle='-', color='blue', label='Terminal Velocity')

# 그래프 설정
plt.title("Terminal Velocity by Drop Radius", fontsize=14)
plt.xlabel("Drop Radius (m)", fontsize=12)
plt.ylabel("Terminal Velocity (m/s)", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.tight_layout()

# 플롯 표시
plt.show()
