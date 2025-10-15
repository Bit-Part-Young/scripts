#!/usr/bin/env python3

"""绘制单条 GSFE 曲线"""

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


def plot_gsfe(
    gsfe_data_fn: str, label: str, smooth: bool = False, figure_fn: str = "gsfe.png"
):
    """绘制单条 GSFE 曲线"""

    set_roman_plot_params(
        lines_linewidth=3.0, legend_frameon=False, legend_handletextpad=0.2
    )

    fig, ax = plt.subplots(figsize=(8, 6))

    gsfe_data = np.loadtxt(gsfe_data_fn, skiprows=1)
    x, y = gsfe_data[:, 0], gsfe_data[:, 1]

    if not smooth:
        ax.plot(x, y, "o-", label=label)
    else:
        ax.scatter(x, y, marker="o", label=label, s=50)

        # plot smooth curve
        smooth_x, smooth_y = smooth_plot(x, y)
        ax.plot(smooth_x, smooth_y)

    ax.set(
        xlim=(-0.02, 1.02),
        xlabel="$u$/b",
        ylabel="$\gamma$ mJ/m$^2$",
    )

    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    y_max = y.max()
    if y_max <= 500.0:
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(50))
    else:
        ax.yaxis.set_major_locator(MultipleLocator(200))
        ax.yaxis.set_minor_locator(MultipleLocator(100))

    y_min = y.min()
    if y_min < 0.0:
        ax.set_ylim(bottom=-30.0)
    else:
        ax.set_ylim(bottom=-20.0)

    ax.legend()

    fig.savefig(figure_fn)

    print(f"\n{figure_fn} is generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot GSFE curve.", epilog="Author: SLY."
    )

    parser.add_argument("gsfe_data_fn", default="gsfe.dat", help="GSFE data filename")
    parser.add_argument("figure_fn", default="gsfe.png", help="output figure filename")
    parser.add_argument("-l", "--label", metavar="STR", help="label")
    parser.add_argument("-s", "--smooth", action="store_true", help="smooth the curve")

    args = parser.parse_args()

    plot_gsfe(args.gsfe_data_fn, args.label, args.smooth, args.figure_fn)
