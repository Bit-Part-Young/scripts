#!/usr/bin/env python3

"""
根据 VTST 生成的 neb.dat spline.dat 文件绘制 NEB 曲线

数据文件生成: nebef.pl; nebbarrier.pl; nebspline.pl
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d

from spt.plot_params import set_plot_params


def plot_neb_vtst():
    """根据 VTST 生成的 neb.dat spline.dat 文件绘制 NEB 曲线"""

    neb_data = np.loadtxt("neb.dat")
    distances = neb_data[:, 1]
    energies = neb_data[:, 2]

    distances_normalized = distances / distances[-1]

    # spline interpolation
    if os.path.exists("spline.dat"):
        data_spline = np.loadtxt("spline.dat")
        distances_spline = data_spline[:, 0]
        distances_spline = distances_spline / distances_spline[-1]
        energies_spline = data_spline[:, 2]
    else:
        spline = interp1d(distances_normalized, energies, kind="cubic")
        distances_spline = np.linspace(0, 1, num=100, endpoint=True)
        energies_spline = spline(distances_spline)

    # plotting
    set_plot_params(roman_params=True, sci_params=True)
    fig, ax = plt.subplots()

    ax.plot(distances_normalized, energies, "o", color="#1f77b4", label="NEB")
    ax.plot(distances_spline, energies_spline, "-", color="#1f77b4")

    ax.set_xlabel("Reaction coordinate (normalized)")
    ax.set_ylabel("Energy (eV)")

    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    ax.grid(color="gainsboro", ls="--", lw=0.6)

    ax.legend()

    fig.savefig("neb_vtst.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Plot NEB curve from VTST generated data.")

    args = parser.parse_args()

    plot_neb_vtst()
