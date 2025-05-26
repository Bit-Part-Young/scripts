"""根据 phonopy-bandplot --gnuplot > band.dat 绘制声子谱（示例代码）"""

import matplotlib.pyplot as plt
import pandas as pd
from spt.plot_params import set_plot_params

data = pd.read_csv("pho.dat", sep=" ", header=None, skiprows=2)

set_plot_params()
fig, ax = plt.subplots()

# 以下方式可避免出现连接原点的线
npoints = 51
for i in range(data.shape[0] // npoints):
    ax.plot(
        data.iloc[i * npoints : (i + 1) * npoints, 0],
        data.iloc[i * npoints : (i + 1) * npoints, 1],
        color="black",
    )

# *.dat 文件第 2 行
xtics = [0.00000000, 0.24746350, 0.33495510, 0.59742970, 0.81173940]
# band.yaml 文件中的 band_labels key
xticklabels = ["$\Gamma$", "X", "U|K", "$\Gamma$", "L"]
ax.set(
    xlabel="Wave vector",
    ylabel="Frequency (THz)",
    xlim=(data.iloc[:, 0].min(), data.iloc[:, 0].max()),
    xticks=xtics,
    xticklabels=xticklabels,
)

# 添加水平线和垂直线
ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)
for x in xtics:
    ax.axvline(x, color="gray", linestyle="--", linewidth=0.5)


fig.savefig("phonon_matplotlib.png")
