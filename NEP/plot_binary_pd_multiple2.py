#!/usr/bin/env python3

"""
绘制多条 0K 二元相图（根据数据文件中的 concentration 和 formation energy 数据列绘制）
数据文件的前 2 列为元素浓度（元素顺序需保持一致），最后 1 列为形成能

用于绘制 无 IMC 形成的二元不同 BCC 衍生结构的 0K 相图 DFT、NEP 对比
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from convex_hull import ConvexHull
from matplotlib.axes import Axes
from matplotlib.ticker import MultipleLocator
from spt.plot_params import set_plot_params


def plot_binary_pd(
    ax: Axes, data_fn: str, columns: list[int], linecolor: str, label: str
):
    """绘制 0K 二元相图"""

    # 不添加 structure label，效果不好
    if len(columns) != 2:
        raise ValueError("columns must be a list of 2 integers!")

    df = pd.read_csv(data_fn, sep=None, engine="python")
    element = df.columns[columns[0]]
    # df 按照 element 的值从小到大排序；排序后 B32 都排在 B2 后面，故不再继续调整
    df = df.sort_values(by=element)
    print(df)

    df_concentrations = df.iloc[:, columns[0]]
    df_energies = df.iloc[:, columns[-1]]

    convex_hull = ConvexHull(df_concentrations, df_energies)

    ax.plot(
        convex_hull.concentrations,
        convex_hull.energies,
        "--",
        color=linecolor,
        linewidth=2.0,
        label=label,
    )

    # marker 中心内容不填充，边缘颜色与线条颜色一致
    # 两端用 o 圆 代表纯元素，次两端用 ^ 三角 代表 D03，B2 用 s 方形，B32 用 D 菱形 表示
    # marker 含义后面自己在图题中说明
    marker_styles = ["o", "^", "s", "D", "^", "o"]
    for i in range(len(df_concentrations)):
        ax.scatter(
            df_concentrations.iloc[i],
            df_energies.iloc[i],
            marker=marker_styles[i],
            facecolors="none",
            edgecolors=linecolor,
            linewidths=2.0,
        )

    ax.set_ylim(ymin=-0.2, ymax=0.2)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_xlabel(f"Composition ({element})")
    ax.set_ylabel(f"Formation Energy (eV/atom)")

    return ax


def plot_multiple(data_fn_list: list[str], label_list: list[str], figure_fn: str):
    """绘制多条 0K 二元相图"""

    linecolor_list = [
        "#1f77b4",  # 蓝色 DFT
        "#d62728",  # 红色 NEP
        "#2ca02c",  # 绿色 EAM
        "#ff7f0e",  # 橙色
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]

    set_plot_params(roman_params=True, sci_params=True)
    fig, ax = plt.subplots(figsize=(6, 6))

    # cols=[0, -1] 的含义为：数据文件的前 2 列为元素浓度（元素顺序需保持一致），最后 1 列为形成能
    # [ ] 是否将第 1 个相关数据的 markersize 加大
    for i, (data_fn, label) in enumerate(zip(data_fn_list, label_list)):
        ax = plot_binary_pd(
            ax=ax,
            data_fn=data_fn,
            columns=[0, -1],
            linecolor=linecolor_list[i],
            label=label,
        )

    ax.legend(loc="upper right")

    plt.savefig(figure_fn)

    print(f"\n0K binary phase diagram {figure_fn} is generated!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot multiple 0K binary phase diagrams from concentration and formation energy data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
