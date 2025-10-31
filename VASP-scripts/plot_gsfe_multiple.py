#!/usr/bin/env python3

"""绘制多条 GSFE 曲线"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import make_interp_spline

from spt.plot_params import set_roman_plot_params


def smooth_plot(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """平滑曲线"""

    smooth_x = np.linspace(x.min(), x.max(), 300)

    spl = make_interp_spline(x, y, k=3)
    smooth_y = spl(smooth_x)

    return smooth_x, smooth_y


def plot_gsfe_multiple(
    gsfe_fn_list: list[str],
    label_list: list[str],
    smooth: bool = False,
    xlabel: str = "$u$/b",
    figure_fn: str = "gsfe_multiple.png",
):
    """绘制多条 GSFE 曲线"""

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

    # matplotlib 默认颜色
    color_list = [
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
    linestyle_list = [
        "o-",  # 圆形
        ">-",  # 又三角
        "s-",  # 正方形
        "D-",  # 菱形
        "p-",
        "x-",
        "v-",
    ]
    marker_list = ["o", ">", "s", "D", "p", "x", "v"]

    y_max_list = []
    y_min_list = []
    for i, (gsfe_fn, label) in enumerate(zip(gsfe_fn_list, label_list)):

        gsfe_data = np.loadtxt(gsfe_fn, skiprows=1)
        x, y = gsfe_data[:, 0], gsfe_data[:, 1]

        y_max_list.append(y.max())
        y_min_list.append(y.min())

        if not smooth:
            ax.plot(x, y, linestyle_list[i], color=color_list[i], label=label)
        else:
            ax.scatter(
                x, y, marker=marker_list[i], color=color_list[i], label=label, s=50
            )

            smooth_x, smooth_y = smooth_plot(x, y)
            ax.plot(smooth_x, smooth_y, color=color_list[i])

    # xlabel 写法示例: Displacement [111]/b; Slip Distance (<111>/2); Displacement-[112]/2
    ax.set(xlim=(-0.02, 1.02), xlabel=xlabel, ylabel="SFE $\gamma$ (mJ/m$^2$)")

    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    y_max = max(y_max_list)
    if y_max <= 500.0:
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(50))
    else:
        ax.yaxis.set_major_locator(MultipleLocator(200))
        ax.yaxis.set_minor_locator(MultipleLocator(100))

    y_min = min(y_min_list)
    if y_min < 0.0:
        ax.set_ylim(bottom=-30.0)
    else:
        ax.set_ylim(bottom=-20.0)

    ax.legend()

    fig.savefig(figure_fn)

    # 关闭图形以释放内存
    plt.close(fig)

    print(f"\n{figure_fn} is generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot multiple GSFE curves.", epilog="Author: SLY."
    )

    parser.add_argument(
        "-i", "--gsfe_fn_list", nargs="+", metavar="FILE", help="GSFE data filenames"
    )
    parser.add_argument("-l", "--label_list", nargs="+", metavar="STR", help="labels")
    parser.add_argument("-xl", "--xlabel", metavar="STR", help="xlabel")
    parser.add_argument(
        "-o",
        "--figure_fn",
        metavar="FILE",
        default="gsfe_multiple.png",
        help="output figure filename",
    )
    parser.add_argument("-s", "--smooth", action="store_true", help="smooth the curve")

    args = parser.parse_args()

    plot_gsfe_multiple(
        args.gsfe_fn_list, args.label_list, args.smooth, args.xlabel, args.figure_fn
    )
