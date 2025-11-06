#!/usr/bin/env python3

"""绘制 0K 二元相图（根据数据文件中的 concentration 和 formation energy 数据列绘制）"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from convex_hull import ConvexHull
from matplotlib.ticker import MultipleLocator

from spt.plot_params import set_plot_params


def plot_binary_pd(data_fn: str, cols: list[int], figure_fn: str):
    """绘制 0K 二元相图（根据数据文件中的 concentration 和 formation energy 数据列绘制）"""

    df = pd.read_csv(data_fn, sep=None, engine="python")
    df_concentrations = df.iloc[:, cols[0]]
    element = df.columns[cols[0]]
    df_energies = df.iloc[:, cols[1]]

    convex_hull = ConvexHull(df_concentrations, df_energies)

    set_plot_params(roman_params=True, sci_params=True)
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(convex_hull.concentrations, convex_hull.energies, color="black")

    ax.scatter(df_concentrations, df_energies)

    # [ ] y 轴范围设置待优化
    ax.set_ylim(top=0.3)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_xlabel(f"Composition ({element})")
    ax.set_ylabel(f"Formation Energy (eV/atom)")

    fig.savefig(figure_fn)

    print(f"\n0K binary phase diagram {figure_fn} is generated!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot binary 0K phase diagram from concentration and formation energy data.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "data_fn",
        help="data filename (must contain concentration and formation energy columns)",
    )
    parser.add_argument(
        "cols",
        type=int,
        nargs=2,
        help="column numbers for concentration and formation energy (starting from 0)",
    )
    parser.add_argument(
        "-o",
        "--fig_fn",
        default="binary_pd.png",
        metavar="FILE",
        help="output figure filename",
    )

    args = parser.parse_args()

    plot_binary_pd(args.data_fn, args.cols, args.fig_fn)
