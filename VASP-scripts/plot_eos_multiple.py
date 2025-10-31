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

    # 使所有的字体加粗
    # 增加线、轴的线宽（3.0 有加粗效果）
    # 增加主、次刻度线长度与宽度
    # 不显示图例边框
    set_roman_plot_params(
        font_weight="bold",
        axes_labelweight="bold",
        axes_linewidth=3.0,
        lines_linewidth=3.0,
        xtick_major_width=3.0,
        ytick_major_width=3.0,
        xtick_minor_width=3.0,
        ytick_minor_width=3.0,
        xtick_major_size=6.0,
        ytick_major_size=6.0,
        xtick_minor_size=3.5,
        ytick_minor_size=3.5,
        legend_frameon=False,
    )

    fig, ax = plt.subplots(figsize=(8, 6))

    # 使得 axes 框为正方形
    # ax.set_box_aspect(1)

    prop_cycle_list = [
        "#1f77b4",  # 蓝色
        "#d62728",  # 红色
        "#2ca02c",  # 绿色
        "#ff7f0e",  # 橙色
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
    linestyle_list = [
        "o-",  # 圆形
        ">-",  # 又三角
        "s-",  # 正方形
        "D-",  # 菱形
        "p-",
        "x-",
        "v-",
    ]

    for i, (data_fn, label) in enumerate(zip(data_fn_list, label_list)):
        df = pd.read_csv(data_fn, sep=None, engine="python")

        if not deltaE:
            df_volume = df.iloc[:, 1]
            df_energy = df.iloc[:, 2]
        else:
            df_volume = df.iloc[:, 0]
            df_energy = df.iloc[:, 1]

        if i == 0:
            ax.plot(
                df_volume,
                df_energy,
                linestyle_list[i],
                markersize=10,
                color=prop_cycle_list[i],
                label=label,
            )
        else:
            ax.plot(
                df_volume,
                df_energy,
                linestyle_list[i],
                markersize=7,
                markerfacecolor="none",
                color=prop_cycle_list[i],
                label=label,
            )

    ax.legend()

    ax.set_xlim(int(df_volume.min()), int(df_volume.max()) + 1.1)

    ax.xaxis.set_minor_locator(MultipleLocator(1.0))

    ax.set_xlabel("Volume ($\AA^3$/atom)")
    if not deltaE:
        ax.set_ylabel("Energy (eV/atom)")
    else:
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.set_ylabel("$\Delta E$ (eV/atom)")

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
