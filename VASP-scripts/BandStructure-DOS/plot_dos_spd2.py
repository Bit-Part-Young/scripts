"""体系分态密度绘制（不用 pymatgen 模块）"""

import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from spt.plot_params import set_roman_plot_params

csv_fname = "dos_spd.csv"
df = pd.read_csv(csv_fname)
energy = df["energy"]

set_roman_plot_params()
fig, ax = plt.subplots()

for i in range(1, df.shape[1]):
    # 平滑处理
    # dos_smooth = gaussian_filter1d(df.iloc[:, i], sigma=0.05)
    # ax.plot(energy, dos_smooth, label=df.columns[i])
    ax.plot(energy, df.iloc[:, i], label=df.columns[i])

ax.axvline(x=0, color="black", linestyle="--")

ax.set_xlim(-10, 2.5)
ax.set_ylim(0, 35)

ax.set_xlabel("Energy (eV)")
ax.set_ylabel("DOS")

ax.legend()

fig.savefig("dos_spd.png")

print("DOS Figure is generated.")
