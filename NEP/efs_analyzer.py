#!/usr/bin/env python3

"""绘制 energy, force, virial/stress, natoms 直方图"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes


def forces_histogram(data_fn: str, property_name: str, bins: int = 30):
    """绘制 force 直方图"""

    array = np.loadtxt(data_fn)
    label_list = ["fx", "fy", "fz"]

    df = pd.DataFrame(array, columns=label_list)
    pd.set_option("display.float_format", "{:.2f}".format)
    print()
    print(df.describe())

    force_max_row = df.max(axis=1).abs()
    count1 = ((force_max_row >= 0.0) & (force_max_row <= 5.0)).sum()
    count2 = (force_max_row > 5.0).sum()
    count3 = (force_max_row > 10.0).sum()
    count4 = (force_max_row > 15.0).sum()
    count5 = (force_max_row > 20.0).sum()

    print()
    print(f"force (abs) in 0.0~5.0: {count1}")
    print(f"force (abs) > 5.0: {count2}")
    print(f"force (abs) > 10.0: {count3}")
    print(f"force (abs) > 15.0: {count4}")
    print(f"force (abs) > 20.0: {count5}")

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 4), dpi=500)

    for i, ax in enumerate(axs.flat):
        ax: Axes
        ax.hist(array[:, i], bins=bins, edgecolor="black", label=label_list[i])
        ax.set_xlabel(f"{label_list[i]}")
        ax.set_ylabel("Frequency")

        ax.legend()

    plt.tight_layout()

    fig_fn = f"range_{property_name}.png"
    fig.savefig(fig_fn)

    print(f"\nFigure {fig_fn} generated!")


def virial_histogram(data_fn: str, property_name: str, bins: int = 30):
    """绘制 virial/stress 直方图"""

    array = np.loadtxt(data_fn)
    label_list = ["xx", "yy", "zz", "xy", "xz", "yz"]

    df = pd.DataFrame(array, columns=label_list)
    pd.set_option("display.float_format", "{:.2f}".format)
    print()
    print(df.describe())
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 12), dpi=500)

    for i, ax in enumerate(axs.flat):
        ax: Axes
        ax.hist(array[:, i], bins=bins, edgecolor="black", label=label_list[i])
        ax.set_xlabel(f"{label_list[i]}")
        ax.set_ylabel("Frequency")

        ax.legend()

    plt.tight_layout()

    fig_fn = f"range_{property_name}.png"
    fig.savefig(fig_fn)

    print(f"\nFigure {fig_fn} generated!")


def energy_histogram(data_fn: str, property_name: str, bins: int = 30):
    """绘制 energy, natoms 直方图"""

    array = np.loadtxt(data_fn)

    df = pd.DataFrame(array, columns=[property_name])
    pd.set_option("display.float_format", "{:.2f}".format)
    print()
    print(df.describe())

    fig, ax = plt.subplots(figsize=(6, 4), dpi=300)

    ax.hist(array, bins=bins, edgecolor="black", label=property_name)

    ax.set_title(f"{property_name} Histogram")
    ax.set_xlabel(f"{property_name}")
    ax.set_ylabel("Frequency")

    ax.legend()

    plt.tight_layout()

    fig_fn = f"range_{property_name}.png"
    fig.savefig(fig_fn)

    print(f"\nFigure {fig_fn} generated!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Histogram plot of energy, forces, virial/stress, natoms.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument("data_fn", type=str, help="data filename")

    parser.add_argument("--bins", type=int, default=30, help="bins number")
    args = parser.parse_args()

    property_name = args.data_fn.split("_")[0]
    if property_name in ["energy", "natoms"]:
        energy_histogram(args.data_fn, property_name, args.bins)
    elif property_name in ["virial", "stress"]:
        virial_histogram(args.data_fn, property_name, args.bins)
    elif property_name == "force":
        forces_histogram(args.data_fn, property_name, args.bins)
    else:
        raise ValueError(f"Invalid property name: {property_name}")
