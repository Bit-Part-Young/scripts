#!/usr/bin/env python3

"""
绘制 NEP 训练能量、力、virial/stress DFT 计算值 与 NEP 预测值对比图
训练集、测试集数据绘制在一起
"""


import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from spt.plot_params import set_plot_params


def load_data(
    energy_train_file: str = "energy_train.out",
    force_train_file: str = "force_train.out",
    virial_train_file: str = "virial_train.out",
    stress_train_file: str = "stress_train.out",
    energy_test_file: str = "energy_test.out",
    force_test_file: str = "force_test.out",
    virial_test_file: str = "virial_test.out",
    stress_test_file: str = "stress_test.out",
):
    """加载数据文件"""

    energy_train = np.loadtxt(energy_train_file)
    force_train = np.loadtxt(force_train_file)
    virial_train = np.loadtxt(virial_train_file)
    stress_train = np.loadtxt(stress_train_file)

    energy_test = np.loadtxt(energy_test_file)
    force_test = np.loadtxt(force_test_file)
    virial_test = np.loadtxt(virial_test_file)
    stress_test = np.loadtxt(stress_test_file)

    return (
        energy_train,
        force_train,
        virial_train,
        stress_train,
        energy_test,
        force_test,
        virial_test,
        stress_test,
    )


def calculate_rmse(
    predicted: np.ndarray,
    dft: np.ndarray,
):
    """计算 RMSE"""

    rmse = np.sqrt(np.mean((predicted - dft) ** 2))

    return rmse


def parity_plot(
    ax: Axes,
    dft_values,
    predicted_values,
    color,
    label,
):
    """绘制 DFT 计算值与 NEP 预测值对比图"""

    ax.scatter(
        dft_values,
        predicted_values,
        s=50,
        c=color,
        label=f"{label}",
    )


def plot_training_test(
    stress: bool = False,
):
    """绘制能量、力和应力的比较图"""

    (
        energy_train,
        force_train,
        virial_train,
        stress_train,
        energy_test,
        force_test,
        virial_test,
        stress_test,
    ) = load_data()

    set_plot_params(
        legend_fontsize=20,
        legend_labelspacing=0.2,
        savefig_dpi=500,
        axes_grid=True,
    )

    plt.figure(figsize=(20, 10))

    text_fontsize = 17

    # ------------------  能量对比图 + 误差直方图  ------------------
    plt.subplot(1, 3, 1)
    ax = plt.gca()

    parity_plot(
        ax,
        energy_train[:, 1],
        energy_train[:, 0],
        "#2472A3",
        "Train",
    )

    parity_plot(
        ax,
        energy_test[:, 1],
        energy_test[:, 0],
        "#EECA40",
        "Test",
    )

    energy_train_rmse = calculate_rmse(energy_train[:, 0], energy_train[:, 1]) * 1000
    ax.text(
        0.3,
        0.15,
        f"Train RMSE: {energy_train_rmse:.2f} meV/$\mathrm{{\AA}}$",
        transform=ax.transAxes,
        fontsize=text_fontsize,
    )

    energy_test_rmse = calculate_rmse(energy_test[:, 0], energy_test[:, 1]) * 1000
    ax.text(
        0.3,
        0.08,
        f"Test RMSE: {energy_test_rmse:.2f} meV/$\mathrm{{\AA}}$",
        transform=ax.transAxes,
        fontsize=text_fontsize,
    )

    xmin, xmax = ax.get_xlim()
    ax.plot([xmin, xmax], [xmin, xmax], c="black", lw=1)

    ax.set(
        xlim=[xmin, xmax],
        ylim=[xmin, xmax],
        xlabel="DFT energy (eV/atom)",
        ylabel="NEP energy (eV/atom)",
    )

    ax.legend(loc="upper left")

    ax.set_aspect("equal", "box")

    # ------------------  力对比图 + 误差直方图  ------------------
    plt.subplot(1, 3, 2)
    ax = plt.gca()

    parity_plot(
        ax,
        force_train[:, 3:].flatten(),
        force_train[:, :3].flatten(),
        "#2472A3",
        "Train",
    )

    parity_plot(
        ax,
        force_test[:, 3:].flatten(),
        force_test[:, :3].flatten(),
        "#EECA40",
        "Test",
    )

    force_train_rmse = calculate_rmse(force_train[:, 0:3], force_train[:, 3:6]) * 1000
    ax.text(
        0.3,
        0.15,
        f"Train RMSE: {force_train_rmse:.2f} meV/$\mathrm{{\AA}}$",
        transform=ax.transAxes,
        fontsize=text_fontsize,
    )

    force_test_rmse = calculate_rmse(force_test[:, 0:3], force_test[:, 3:6]) * 1000
    ax.text(
        0.3,
        0.08,
        f"Test RMSE: {force_test_rmse:.2f} meV/$\mathrm{{\AA}}$",
        transform=ax.transAxes,
        fontsize=text_fontsize,
    )

    xmin, xmax = ax.get_xlim()
    ax.plot([xmin, xmax], [xmin, xmax], c="black", lw=1)

    ax.set(
        xlim=[xmin, xmax],
        ylim=[xmin, xmax],
        xlabel="DFT force (eV/$\mathrm{\AA}$)",
        ylabel="NEP force (eV/$\mathrm{\AA}$)",
    )

    ax.legend(loc="upper left")

    ax.set_aspect("equal", "box")

    # ------------------  位力/应力对比图 + 误差直方图  ------------------
    plt.subplot(1, 3, 3)
    ax = plt.gca()

    if stress:
        virial_train = stress_train
        virial_test = stress_test
        axes_label = "stress (GPa)"
    else:
        axes_label = "virial (eV/atom)"

    parity_plot(
        ax,
        virial_train[:, 6:].flatten(),
        virial_train[:, :6].flatten(),
        "#2472A3",
        "Train",
    )

    parity_plot(
        ax,
        virial_test[:, 6:].flatten(),
        virial_test[:, :6].flatten(),
        "#EECA40",
        "Test",
    )

    virial_train_rmse = calculate_rmse(virial_train[:, 0:6], virial_train[:, 6:]) * 1000
    ax.text(
        0.3,
        0.15,
        f"Train RMSE: {virial_train_rmse:.2f} meV/$\mathrm{{\AA}}$",
        transform=ax.transAxes,
        fontsize=text_fontsize,
    )

    virial_test_rmse = calculate_rmse(virial_test[:, 0:6], virial_test[:, 6:]) * 1000
    ax.text(
        0.3,
        0.08,
        f"Test RMSE: {virial_test_rmse:.2f} meV/$\mathrm{{\AA}}$",
        transform=ax.transAxes,
        fontsize=text_fontsize,
    )

    xmin, xmax = ax.get_xlim()
    ax.plot([xmin, xmax], [xmin, xmax], c="black", lw=1)

    ax.set(
        xlim=[xmin, xmax],
        ylim=[xmin, xmax],
        xlabel=f"DFT {axes_label}",
        ylabel=f"NEP {axes_label}",
    )

    ax.legend(loc="upper left")

    ax.set_aspect("equal", "box")

    # ------------------  保存图片  ------------------
    plt.tight_layout()
    plt.savefig("nep_train_test.png", bbox_inches="tight")

    print("\nNEP train & test results figure is generated!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot NEP train & test results (parity plot of energy, force, virial, stress).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--stress",
        action="store_true",
        help="parity plot of stress or virial",
    )

    args = parser.parse_args()

    plot_training_test(args.stress)
