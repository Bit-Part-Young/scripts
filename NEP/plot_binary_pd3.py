#!/usr/bin/env python3

"""
绘制 0K 二元相图（根据数据文件中的 concentration, structure label, formation energy 数据列绘制）
添加 结构 label 文本到每个数据点上（文本可能会与数据点有重叠）

用于绘制 无 IMC 形成的二元不同 BCC 衍生结构的 0K 相图
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from convex_hull import ConvexHull
from matplotlib.ticker import MultipleLocator

from spt.plot_params import set_plot_params


def plot_binary_pd(data_fn: str, columns: list[int], figure_fn: str):
    """绘制 0K 二元相图"""

    if len(columns) != 3:
        raise ValueError("columns must be a list of 3 integers!")

    df = pd.read_csv(data_fn, sep=None, engine="python")
    element = df.columns[columns[0]]
    df_concentrations = df.iloc[:, columns[0]]
    df_energies = df.iloc[:, columns[2]]
    # 添加 structure 标签 如 D03 B32 B2 GS 等
    df_structure_labels = df.iloc[:, columns[1]]

    convex_hull = ConvexHull(df_concentrations, df_energies)

    color = "#1f77b4"  # 蓝色

    set_plot_params(roman_params=True, sci_params=True)
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(convex_hull.concentrations, convex_hull.energies, "o--", color=color)

    ax.scatter(df_concentrations, df_energies, marker="s")

    ax.set_ylim(ymin=-0.2, ymax=0.2)

    # 在每个点上添加 structure label
    # 计算合适的垂直偏移量（基于 y 轴范围）
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    # 偏移量为y轴范围的 2%
    text_offset = y_range * 0.02
    for i in range(len(df_concentrations)):
        # 若 concentration 接近于 0 或 1，则不添加 label
        if df_concentrations[i] < 0.01 or df_concentrations[i] > 0.99:
            continue

        ax.text(
            df_concentrations[i],
            df_energies[i] + text_offset,
            df_structure_labels[i],
            fontsize=14,
            color=color,
            # ha="center",  # 水平居中对齐
            va="bottom",  # 垂直底部对齐，配合偏移量使文字在点上方
        )

    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_xlabel(f"Composition ({element})")
    ax.set_ylabel(f"Formation Energy (eV/atom)")

    fig.savefig(figure_fn)

    print(f"\n0K binary phase diagram {figure_fn} is generated!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot binary 0K phase diagram from concentration, structure label and formation energy data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument("data_fn", help="data filename")
    parser.add_argument(
        "-c",
        "--columns",
        type=int,
        nargs=3,
        metavar="N",
        help="column numbers for concentration, structure label and formation energy (starting from 0)",
    )
    parser.add_argument(
        "-o",
        "--fig_fn",
        default="binary_pd.png",
        metavar="FILE",
        help="output figure filename",
    )

    args = parser.parse_args()

    plot_binary_pd(args.data_fn, args.columns, args.fig_fn)
