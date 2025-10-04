#!/usr/bin/env python3

"""EOS 曲线多条绘制（支持 e-v 和 deltaE-v 数据）"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MultipleLocator
from spt.plot_params import set_roman_plot_params


def plot_eos_multiple(
    data_fn_list: list[str], label_list: list[str], figure_fn: str, deltaE: bool = False
):
    """EOS 曲线多条绘制"""

    set_roman_plot_params(lines_linewidth=3.0)
    fig, ax = plt.subplots(figsize=(8, 6))

    linestyle_list = ["o-", "s-", "D-", "p-", "x-", "v-"]

    for i, (data_fn, label) in enumerate(zip(data_fn_list, label_list)):
        df = pd.read_csv(data_fn, sep=None, engine="python")

        if not deltaE:
            df_volume = df.iloc[:, 1]
            df_energy = df.iloc[:, 2]
        else:
            df_volume = df.iloc[:, 0]
            df_energy = df.iloc[:, 1]

        ax.plot(
            df_volume, df_energy, linestyle_list[i], markerfacecolor="none", label=label
        )

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
        description="Plot multiple EOS curves.", epilog="Author: SLY."
    )

    parser.add_argument(
        "-i",
        "--data_fn_list",
        nargs="+",
        metavar="FILE",
        help="input EOS data filenames",
    )
    parser.add_argument("-l", "--label_list", nargs="+", metavar="STR", help="labels")
    parser.add_argument(
        "-o",
        "--figure_fn",
        default="eos.png",
        metavar="FILE",
        help="output figure filename",
    )
    parser.add_argument(
        "--deltaE", action="store_true", help="whether plot deltaE-v EOS curve"
    )

    args = parser.parse_args()

    plot_eos_multiple(args.data_fn_list, args.label_list, args.figure_fn, args.deltaE)
