#!/usr/bin/env python3

"""EOS 曲线绘制（支持 e-v 和 deltaE-v 数据）"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MultipleLocator

from spt.plot_params import set_roman_plot_params


def plot_eos(
    data_fn: str = "eos.dat", figure_fn: str = "eos.png", deltaE: bool = False
):
    """EOS 曲线绘制"""

    df = pd.read_csv(data_fn, sep=None, engine="python")

    if not deltaE:
        df_volume = df.iloc[:, 1]
        df_energy = df.iloc[:, 2]
    else:
        df_volume = df.iloc[:, 0]
        df_energy = df.iloc[:, 1]

    set_roman_plot_params(lines_linewidth=3.0)
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(df_volume, df_energy, "o-", markerfacecolor="none", label="EOS")

    ax.legend()

    ax.set_xlim(int(df_volume.min()), int(df_volume.max()) + 1.1)

    ax.xaxis.set_minor_locator(MultipleLocator(1.0))

    ax.set_xlabel("Volume ($\AA^3$/atom)")
    if not deltaE:
        ax.set_ylabel("Energy (eV/atom)")
    else:
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.set_ylabel(f"$\Delta$E (eV/atom)")

    fig.savefig(figure_fn)

    print(f"\n{figure_fn} is generated.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot EOS curve.", epilog="Author: SLY."
    )

    parser.add_argument(
        "data_fn", nargs="?", default="eos.dat", help="input EOS data filename"
    )
    parser.add_argument("figure_fn", default="eos.png", help="output figure filename")
    parser.add_argument(
        "--deltaE", action="store_true", help="whether plot deltaE-v EOS curve"
    )

    args = parser.parse_args()

    plot_eos(args.data_fn, args.figure_fn, args.deltaE)
