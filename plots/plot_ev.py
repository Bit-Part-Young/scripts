"""
EV（能量-体积）曲线绘制
数据文件内容:
晶格常数、平均原子体积、平均原子能量
a volume_pa energy_pa
"""

import matplotlib.pyplot as plt
import pandas as pd
from spt.plot_params import set_roman_plot_params


def plot_ev(
    csv_fn: str,
    output_fn: str = "ev.png",
):
    """EV（能量-体积）曲线绘制"""

    df = pd.read_csv(
        csv_fn,
        sep=None,
        engine="python",
    )

    df_volume = df.iloc[:, 1]
    df_energy = df.iloc[:, 2]

    set_roman_plot_params()
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(
        df_volume,
        df_energy,
        marker="o",
        markerfacecolor="none",
        linestyle="-",
    )

    ax.set(
        xlabel="Volume ($\AA^3$/atom)",
        ylabel="Energy (eV/atom)",
    )

    fig.savefig(output_fn)


if __name__ == "__main__":
    csv_fn = "ev.dat"

    plot_ev(csv_fn)
