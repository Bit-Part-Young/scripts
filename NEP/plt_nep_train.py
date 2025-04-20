#!/usr/bin/env python3

"""
绘制 NEP 训练 loss 演化与能量、力、virial/stress DFT 计算值 与 NEP 预测值对比图
训练集、测试集数据分别绘图

reference:
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from spt.plot_params import set_plot_params


def load_data(test: bool = False):
    """加载数据文件"""

    if test:
        label = "test"
    else:
        label = "train"

    loss_data = np.loadtxt(f"loss.out")
    energy_data = np.loadtxt(f"energy_{label}.out")
    force_data = np.loadtxt(f"force_{label}.out")
    virial_data = np.loadtxt(f"virial_{label}.out")
    stress_data = np.loadtxt(f"stress_{label}.out")

    return loss_data, energy_data, force_data, virial_data, stress_data


def calculate_rmse(
    predicted_values: np.ndarray,
    dft_values: np.ndarray,
):
    """计算 RMSE"""

    return np.sqrt(np.mean((predicted_values - dft_values) ** 2))


def calculate_limits(
    data: np.ndarray,
    padding: float = 0.08,
):
    """计算轴动态范围"""

    data_min = np.min(data)
    data_max = np.max(data)
    data_range = data_max - data_min

    return data_min - padding * data_range, data_max + padding * data_range


def parity_plot(
    ax: Axes,
    dft_values: np.ndarray,
    predicted_values: np.ndarray,
    xmin: float,
    xmax: float,
):
    """绘制 DFT 计算值与 NEP 预测值对比图"""

    if dft_values.ndim == 1:
        ax.scatter(dft_values, predicted_values, s=50)
    else:
        for i in range(0, dft_values.shape[1]):
            ax.scatter(dft_values[:, i], predicted_values[:, i], s=50)

    ax.plot([xmin, xmax], [xmin, xmax], color="grey", linestyle="--")


def plot_nep_train(
    test: bool = False,
    stress: bool = False,
):
    """绘制 NEP 训练结果"""

    loss, energy_data, force_data, virial_data, stress_data = load_data(test)

    set_plot_params(
        legend_fontsize=20,
        legend_handletextpad=0.2,
        savefig_dpi=500,
        axes_grid=True,
    )

    fig, axs = plt.subplots(2, 2, figsize=(18, 15))

    # ------------------  损失演化曲线  ------------------
    ax: Axes = axs[0, 0]
    ax.grid(False)

    loss[:, 0] = np.arange(1, len(loss) + 1) * 100
    ax.loglog(loss[:, 0], loss[:, 1:7], "-")

    ax.set(
        xlim=(1e2, loss[:, 0].max()),
        xlabel="Generation",
        ylabel="Loss",
    )

    ax.legend(
        labels=["Total", "L1", "L2", "E-train", "F-train", "V-train"],
        ncols=3,
        loc="lower left",
    )

    # ------------------  能量对比图 + 误差直方图  ------------------
    energy_min, energy_max = calculate_limits(energy_data[:, 1])
    ax: Axes = axs[0, 1]

    parity_plot(
        ax,
        energy_data[:, 1],
        energy_data[:, 0],
        energy_min,
        energy_max,
    )

    ax.set(
        xlim=[energy_min, energy_max],
        ylim=[energy_min, energy_max],
        xlabel="DFT energy (eV/atom)",
        ylabel="NEP energy (eV/atom)",
    )

    ax.legend(labels=["energy"])

    energy_rmse = calculate_rmse(energy_data[:, 0], energy_data[:, 1]) * 1000
    ax.text(
        0.5,
        0.08,
        f"RMSE: {energy_rmse:.2f} meV/atom",
        transform=ax.transAxes,
        fontsize=20,
    )

    # ------------------  力对比图 + 误差直方图  ------------------
    force_min, force_max = calculate_limits(force_data[:, 3:6].reshape(-1))
    ax: Axes = axs[1, 0]

    parity_plot(
        ax,
        force_data[:, 3:6],
        force_data[:, 0:3],
        force_min,
        force_max,
    )

    ax.set(
        xlim=[force_min, force_max],
        ylim=[force_min, force_max],
        xlabel=r"DFT force (eV/$\AA$)",
        ylabel=r"NEP force (eV/$\AA$)",
    )

    ax.legend(labels=["fx", "fy", "fz"], ncols=3)

    force_rmse = calculate_rmse(force_data[:, 0:3], force_data[:, 3:6]) * 1000
    ax.text(
        0.5,
        0.08,
        rf"RMSE: {force_rmse:.2f} meV/$\AA$",
        transform=ax.transAxes,
        fontsize=20,
    )

    # ------------------  位力/应力对比图 + 误差直方图  ------------------
    if stress:
        virial_data = stress_data
        axes_label = "stress (GPa)"
        rmse_unit = "GPa"
    else:
        axes_label = "virial (eV/atom)"
        rmse_unit = "eV/atom"

    virial_min, virial_max = calculate_limits(virial_data[:, 6:12].reshape(-1))
    ax: Axes = axs[1, 1]

    parity_plot(
        ax,
        virial_data[:, 6:12],
        virial_data[:, 0:6],
        virial_min,
        virial_max,
    )

    ax.set(
        xlim=[virial_min, virial_max],
        ylim=[virial_min, virial_max],
        xlabel=f"DFT {axes_label}",
        ylabel=f"NEP {axes_label}",
    )

    ax.legend(
        labels=["xx", "yy", "zz", "xy", "yz", "zx"],
        ncols=3,
    )

    virial_rmse = calculate_rmse(virial_data[:, 0:6], virial_data[:, 6:12])
    ax.text(
        0.5,
        0.08,
        f"RMSE: {virial_rmse:.4f} {rmse_unit}",
        transform=ax.transAxes,
        fontsize=20,
    )

    plt.tight_layout()

    if test:
        fig.savefig("nep_test.png")

        print("\nnep_test.png generated!")
    else:
        fig.savefig("nep_train.png")

        print("\nnep_train.png generated!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot NEP training results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--test",
        action="store_true",
        help="plot training or test sets",
    )

    parser.add_argument(
        "--stress",
        action="store_true",
        help="parity plot of stress or virial",
    )
    args = parser.parse_args()

    plot_nep_train(args.test, args.stress)
