#!/usr/bin/env python3

"""
绘制多条 0K 二元相图（根据数据文件中的 concentration 和 formation energy 数据列绘制）
数据文件的前 2 列为元素浓度（元素顺序需保持一致），最后 1 列为形成能
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from convex_hull import ConvexHull
from matplotlib.axes import Axes
from matplotlib.ticker import MultipleLocator

from spt.plot_params import set_plot_params


def plot_binary_pd(
    ax: Axes, data_fn: str, cols: list[int], linecolor: str, marker: str, label: str
):
    """绘制 0K 二元相图（根据数据文件中的 concentration 和 formation energy 数据列绘制）"""

    df = pd.read_csv(data_fn, sep=None, engine="python")
    df_concentrations = df.iloc[:, cols[0]]
    element = df.columns[cols[0]]
    df_energies = df.iloc[:, cols[1]]

    convex_hull = ConvexHull(df_concentrations, df_energies)

    ax.plot(convex_hull.concentrations, convex_hull.energies, color=linecolor)

    ax.scatter(df_concentrations, df_energies, label=label, marker=marker)

    ax.set_ylim(top=0.3)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_xlabel(f"Composition ({element})")
    ax.set_ylabel(f"Formation Energy (eV/atom)")

    return ax


def plot_multiple(data_fn_list: list[str], label_list: list[str], figure_fn: str):
    """绘制多条 0K 二元相图"""

    # [ ] 优化 marker 样式和线条颜色
    marker_list = ["o", "x", "s", "D", "p", ">", "v"]
    linecolor_list = ["black", "red"]

    set_plot_params(roman_params=True, sci_params=True)
    fig, ax = plt.subplots(figsize=(6, 6))

    # cols=[0, -1] 的含义为：数据文件的前 2 列为元素浓度（元素顺序需保持一致），最后 1 列为形成能
    for i, (data_fn, label) in enumerate(zip(data_fn_list, label_list)):
        ax = plot_binary_pd(
            ax=ax,
            data_fn=data_fn,
            cols=[0, -1],
            linecolor=linecolor_list[i],
            marker=marker_list[i],
            label=label,
        )

    ax.legend()

    plt.savefig(figure_fn)

    print(f"\n0K binary phase diagram {figure_fn} is generated!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot multiple 0K binary phase diagrams from concentration and formation energy data",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "-i", "--data_fn_list", nargs="+", metavar="FILE", help="data filenames"
    )
    parser.add_argument("-l", "--label_list", nargs="+", metavar="STR", help="labels")
    parser.add_argument(
        "-o",
        "--fig_fn",
        default="binary_pd.png",
        metavar="FILE",
        help="output figure filename",
    )

    args = parser.parse_args()

    plot_multiple(args.data_fn_list, args.label_list, args.fig_fn)
