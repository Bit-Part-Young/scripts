"""
绘图配色

reference: https://zhuanlan.zhihu.com/p/488125051
"""

import matplotlib.pyplot as plt
import numpy as np

from spt.plot_params import set_roman_plot_params

colors = [
    "#518CD8",
    "#FD6D5A",
    "#6DC354",
    "#FEB40B",
    "#994487",
    "#443295",
]

x = np.linspace(1, 2 * np.e, 20)

set_roman_plot_params(legend_fontsize=15)
fig, ax = plt.subplots()

for i in range(len(colors)):
    y = np.log(x) + i * 0.2

    ax.plot(x, y, "o-", color=colors[i], lw=1, label=colors[i])

ax.legend()

fig.savefig("color_template.png")
