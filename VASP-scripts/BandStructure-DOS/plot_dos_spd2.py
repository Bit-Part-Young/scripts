"""体系分态密度绘制（不用 pymatgen 模块）"""

import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from spt.plot_params import set_roman_plot_params


def plot_dos_spd(
    filter: bool,
    sigma: float,
    figname: str,
):

    csv_fname = "dos_spd.csv"
    df = pd.read_csv(csv_fname)
    energy = df["energy"]

    set_roman_plot_params()
    fig, ax = plt.subplots()

    for i in range(1, df.shape[1]):
        if filter:
            # 平滑处理
            dos_smooth = gaussian_filter1d(df.iloc[:, i], sigma=sigma)
            ax.plot(
                energy,
                dos_smooth,
                label=df.columns[i],
            )
        else:
            ax.plot(
                energy,
                df.iloc[:, i],
                label=df.columns[i],
            )

    ax.axvline(
        x=0,
        color="black",
        linestyle="--",
    )

    ax.set(
        xlim=(-10, 2.5),
        ylim=(0, 35),
        xlabel="Energy (eV)",
        ylabel="DOS",
    )

    ax.legend(handlelength=1.0)

    fig.savefig(figname)

    print("DOS Figure is generated.")


if __name__ == "__main__":
    plot_dos_spd(
        filter=False,
        figname="dos_spd_raw.png",
    )

    plot_dos_spd(
        filter=True,
        sigma=1.0,
        figname="dos_spd.png",
    )
