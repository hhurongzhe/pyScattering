import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.font_manager import FontProperties


plt.rcParams["font.sans-serif"] = ["PingFang SC"]
# 解决负号'-'显示为方块的问题
plt.rcParams["axes.unicode_minus"] = False

# 设置矩阵尺寸和 walker 数量
matrix_size = 100
num_walkers = 1000000

# 为了结果可复现，设置随机种子
np.random.seed(42)

# 创建一个 100x100 的坐标网格
x, y = np.meshgrid(np.arange(matrix_size), np.arange(matrix_size))

# 创建一个有显著特征的 H0 矩阵作为我们的目标概率分布
# 这里我们创建两个高斯“热点”区域，代表概率较高的地方
peak1 = 2 * np.exp(-((x - 25) ** 2 + (y - 25) ** 2) / (2 * 10**2))
peak2 = 1 * np.exp(-((x - 70) ** 2 + (y - 60) ** 2) / (2 * 15**2))
H0 = peak1 + peak2 + 0.01  # 加入少量背景噪声

# --- 2. 使用“采样-重要性-重采样 (SIR)”生成 Walker ---

print("开始使用 SIR 方法生成 walkers...")

# 步骤 1: 采样 (从简单的均匀分布中生成候选 walker)
# 在 100x100 网格上随机选择 10000 个候选位置
candidate_i = np.random.randint(0, matrix_size, num_walkers)
candidate_j = np.random.randint(0, matrix_size, num_walkers)

# 步骤 2: 加权 (根据 H0 计算每个候选者的重要性权重)
# 权重就是 H0 在候选者位置的值。值越大的地方，权重越高。
weights = H0[candidate_i, candidate_j]

# 步骤 3: 重采样 (根据权重从候选者中抽取最终的 walker)
# 首先，将权重归一化，使其成为一个概率分布
normalized_weights = weights / np.sum(weights)

# 然后，根据归一化后的权重，从候选者中进行有放回的抽取
# 我们抽取其索引，而不是直接抽取坐标
chosen_indices = np.random.choice(np.arange(num_walkers), size=num_walkers, p=normalized_weights, replace=True)

# 得到最终的 walker 坐标
walkers_sir_i = candidate_i[chosen_indices]
walkers_sir_j = candidate_j[chosen_indices]

print("SIR 方法完成！")

# --- 3. (对比方法) 直接从 H0 分布中采样 ---
# 这是在 NumPy 中更直接、高效的方法，可作为验证

print("开始使用直接采样方法生成 walkers...")
# 将 H0 展平为一维数组，并归一化为概率
prob_flat = (H0 / np.sum(H0)).flatten()
# 生成一维索引
indices_flat = np.random.choice(np.arange(matrix_size * matrix_size), size=num_walkers, p=prob_flat, replace=True)
# 将一维索引转换回二维坐标
walkers_direct_i, walkers_direct_j = np.unravel_index(indices_flat, H0.shape)
print("直接采样方法完成！")


# --- 4. 可视化结果 ---

sns.set_theme(style="white")
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
fig.suptitle("重要性采样 (SIR) vs. 直接采样", fontsize=20, y=1.02)


# 图1: 原始目标分布 H0
im = axes[0].imshow(H0, cmap="inferno", origin="lower")
axes[0].set_title("目标分布 H0", fontsize=14)
axes[0].set_xlabel("列 (j)")
axes[0].set_ylabel("行 (i)")

# 图2: SIR 方法生成的 Walker 分布
axes[1].hist2d(walkers_sir_j, walkers_sir_i, bins=matrix_size, cmap="inferno")
axes[1].set_title("重要性采样 (SIR) 生成的分布", fontsize=14)
axes[1].set_xlabel("列 (j)")

# 图3: 直接采样生成的 Walker 分布
axes[2].hist2d(walkers_direct_j, walkers_direct_i, bins=matrix_size, cmap="inferno")
axes[2].set_title("直接采样生成的分布 (验证)", fontsize=14)
axes[2].set_xlabel("列 (j)")

# 为图像添加颜色条
fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6)
plt.tight_layout()
plt.show()
