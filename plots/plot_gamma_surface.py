"""绘制 GSFE 的 γ-surface 图"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from spt.plot_params import set_plot_params

ncols, nrows = 30, 30
df = pd.read_csv("gamma_surface.csv")
X = df.iloc[:, 0].values.reshape(ncols, nrows)
Y = df.iloc[:, 1].values.reshape(ncols, nrows)
Z = df.iloc[:, 2].values.reshape(ncols, nrows)
# eV/Angstrom^2 -> mJ/m^2
Z = Z * 1.6021766208 * 10**4

"""
# 手动生成 GSFE 数据
x = np.linspace(0.0, 1.0, 11)
y = np.linspace(0.0, 1.0, 16)
X, Y = np.meshgrid(x, y)
Z = 1500 * (
    np.sin(np.pi * X) * np.sin(np.pi * Y)
    + 0.5 * np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y)
    + 0.3
)
Z[Z < 0] = 0  # 确保能量非负
"""

set_plot_params(roman_params=True, sci_params=True)
fig, ax = plt.subplots()

# SFE 范围和分段
levels = np.linspace(0, int(Z.max()), 10)

# 填充等高线图
countourf = ax.contourf(X, Y, Z, levels=levels, cmap=plt.get_cmap("jet"))

# 等高线（可选）
contour = ax.contour(X, Y, Z, levels=levels, colors="k", linewidths=1.0, alpha=0.8)

# colorbar
cbar = fig.colorbar(countourf, ax=ax, ticks=levels, label="SFE (mJ/m$^2$)")

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel(r"1/2[$\bar{1}10$]")
ax.set_ylabel(r"1/2[11$\bar{2}$]")
ax.grid(True, linestyle="--", alpha=0.5, which="both")
# ax.set_aspect("equal", adjustable="box")

fig.savefig("gsfe_gamma_surface.png")
