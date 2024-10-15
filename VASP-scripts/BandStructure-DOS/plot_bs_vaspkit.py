"""通过 vaspkit 获取的能带数据进行能带绘制"""

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from spt.plot_params import set_roman_plot_params

band_fn = "REFORMATTED_BAND.dat"
band_df = pd.read_csv(band_fn, sep="\s+", skiprows=1, header=None)

klable_fn = "KLABELS"

klabel_list = []
with open(klable_fn, "r") as f:
    lines = f.readlines()

for line in lines[1:]:
    line = line.strip().split()
    if line:
        klabel_list.append(line)
    else:
        break

xlabel_list = [
    (
        "$\Gamma$"
        if klabel_list[k][0] in ["Gamma", "GAMMA"]
        else klabel_list[k][0]
    )
    for k in range(len(klabel_list))
]

set_roman_plot_params()
fig, ax = plt.subplots()
ax: Axes

for i in range(band_df.shape[1] - 1):
    ax.plot(band_df.iloc[:, 0], band_df.iloc[:, i + 1])

ax.axhline(y=0, color="black", linewidth=0.5)

x_list = []
for j in range(len(klabel_list)):
    x = float(klabel_list[j][1])
    x_list.append(x)
    ax.axvline(x=x, color="black", linestyle="--", linewidth=0.5)

ax.set_xticks(x_list, labels=xlabel_list)

ax.set_xlim(x_list[0], x_list[-1])
ax.set_ylim(-3, 3)

ax.set_ylabel("Energy (eV)")

fig.savefig("bs_vaspkit.jpg")

print("Band structure figure is generated.")
