"""
绘图配色

reference: https://github.com/yh-phys/scripts/blob/master/scripts/python/color.py
"""

import matplotlib.pyplot as plt
import numpy as np

from spt.plot_params import set_roman_plot_params

colors = [
    "#3953a4",
    "#faa316",
    "#d93b2b",
    "#0db14b",
    "#f0a3ff",
    "#0075dc",
    "#993f00",
    "#4c005c",
    "#426600",
    "#ff0010",
    "#9dcc00",
    "#c20088",
    "#003380",
    "#ffa405",
    "#ffff00",
    "#ff5005",
    "#5ef1f2",
    "#740aff",
    "#990000",
    "#00998f",
    "#005c31",
    "#2bce48",
    "#ffcc99",
    "#94ffb5",
    "#8f7c00",
    "#6fa8bb",
    "#808080",
    "#4f4ffe",  # reference: https://zhuanlan.zhihu.com/p/698897001
    "#ce3d32",  # 同上
]

x = np.linspace(0, 2 * np.pi, 700)
y = np.sin(x) + np.cos(x) * np.sin(x) + np.cos(2 * x) ** 2

set_roman_plot_params(
    legend_fontsize=15,
    savefig_dpi=300,
)

fig = plt.figure(figsize=(25, 40))

for i in range(len(colors)):
    ax = fig.add_subplot(len(colors) // 4 + 1, 4, i + 1)
    ax.plot(x, y, color=colors[i], lw=1, label=colors[i])  # alpha=0.9
    ax.fill_between(x, 0, y, facecolor=colors[i])  #  alpha=0.25
    ax.set_xlim(0, 2 * np.pi)
    ax.legend()

fig.savefig("color_template2.png")
