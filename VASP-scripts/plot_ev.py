#!/usr/bin/env python3

"""
EV（能量-体积）曲线绘制
数据文件内容:
晶格常数、平均原子体积、平均原子能量
a volume_pa energy_pa
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd

from spt.plot_params import set_roman_plot_params


def plot_ev(
    data_fn: str = "ev.dat",
    output_fn: str = "ev.png",
):
    """EV（能量-体积）曲线绘制"""

    df = pd.read_csv(data_fn, sep=None, engine="python")

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

    parser = argparse.ArgumentParser(
        description="Plot Energy-Volume curve.",
        epilog="Author: SLY",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "data_fn",
        type=str,
        nargs="?",
        default="ev.dat",
        help="File containing volume (Col2) and energy (Col3) data",
    )

    parser.add_argument(
        "-o",
        "--output_fn",
        type=str,
        nargs="?",
        default="ev.png",
        help="Output figure filename",
    )

    args = parser.parse_args()

    data_fn = args.data_fn
    output_fn = args.output_fn

    plot_ev(data_fn=data_fn, output_fn=output_fn)

    print("Energy-Volume curve figure is generated.")
