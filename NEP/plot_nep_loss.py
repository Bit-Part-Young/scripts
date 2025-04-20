#!/usr/bin/env python3

"""
绘制 NEP 训练 loss 演化与能量、力、virial/stress DFT 计算值 与 NEP 预测值对比图及其对应误差直方图
训练集、测试集数据绘制在一起

reference:
https://github.com/wangchr1617/learning/blob/main/scripts/plot/plot_loss_nep.py
https://github.com/wangchr1617/NEP_GT/blob/main/NEP_Loss/NEP_Loss.ipynb
"""


import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from spt.plot_params import set_plot_params

set_plot_params(
    legend_fontsize=20,
    legend_labelspacing=0.2,
    savefig_dpi=500,
    axes_grid=True,
)


def load_data(
    loss_file: str = "loss.out",
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

    try:
        loss = np.loadtxt(loss_file)
        # 每 100 输出一次；不直接使用 loss[:, 0] 是因为 NEP 训练可能是从续算开始，并非从头开始
        loss[:, 0] = np.arange(1, len(loss) + 1) * 100
        print(f"\nWe have run {loss[-1, 0]:.0f} steps!")
    except Exception as e:
        raise FileNotFoundError(f"\nError loading {loss_file}: {e}")

    # 有时可能没有测试集数据
    def load_optional_file(file_path):
        if os.path.exists(file_path):
            return np.loadtxt(file_path)
        else:
            print(f"Warning: {file_path} not found.")
            return None

    energy_train = np.loadtxt(energy_train_file)
    force_train = np.loadtxt(force_train_file)
    virial_train = np.loadtxt(virial_train_file)
    stress_train = np.loadtxt(stress_train_file)

    energy_test = load_optional_file(energy_test_file)
    force_test = load_optional_file(force_test_file)
    virial_test = load_optional_file(virial_test_file)
    stress_test = load_optional_file(stress_test_file)

    return (
        loss,
        energy_train,
        force_train,
        virial_train,
        stress_train,
        energy_test,
        force_test,
        virial_test,
        stress_test,
    )


def parity_plot(
    ax: Axes,
    dft_values,
    predicted_values,
    color,
    label,
):
    """绘制 DFT 计算值与 NEP 预测值对比图"""

    errors = predicted_values - dft_values
    rmse = np.sqrt(np.mean(errors**2))

    ax.scatter(
        dft_values,
        predicted_values,
        s=50,
        c=color,
        label=f"{label} = {rmse:.4f}",
    )


def plot_hist(
    ax: Axes,
    dft_values,
    predicted_values,
    title,
    bins=15,
    alpha=0.75,
):
    """绘制误差绝对值分布图（归一化）"""

    errors = predicted_values - dft_values
    abs_errors = np.abs(errors)

    inset_ax = ax.inset_axes([0.6, 0.1, 0.35, 0.35])

    inset_ax.hist(
        abs_errors,
        bins=bins,
        color="orange",
        alpha=alpha,
        density=True,
        edgecolor="black",
    )

    inset_ax.set_title(title, fontsize=18)
    inset_ax.set_ylabel("Frequency (%)", fontsize=18)
    inset_ax.set_xlim([0, None])
    inset_ax.set_ylim([0, None])
    inset_ax.set_yticks([])
    inset_ax.yaxis.set_label_position("left")
    inset_ax.yaxis.tick_left()


def plot_training_results(
    stress: bool = False,
    hist: bool = False,
):
    """绘制损失、能量、力和应力的比较图"""

    (
        loss,
        energy_train,
        force_train,
        virial_train,
        stress_train,
        energy_test,
        force_test,
        virial_test,
        stress_test,
    ) = load_data()

    plt.figure(figsize=(18, 15))

    # ------------------  损失演化曲线  ------------------
    plt.subplot(2, 2, 1)
    ax = plt.gca()

    ax.grid(visible=False)

    col_list = [
        "Total",
        r"$L_{1}$",
        r"$L_{2}$",
        "E-train",
        "F-train",
        "V-train",
        "E-test",
        "F-test",
        "V-test",
    ]

    if energy_test is None:
        for i, label in zip(range(1, 7), col_list[:6]):
            ax.loglog(loss[:, 0], loss[:, i], ls="-", lw=2, label=label)
    else:
        for i, label in zip(range(1, 10), col_list):
            ax.loglog(loss[:, 0], loss[:, i], ls="-", lw=2, label=label)

    ax.set(
        xlim=(1e2, loss[:, 0].max()),
        # ylim=(5e-4, 2e0),
        xlabel="Generation",
        ylabel="Loss",
        title="(a) Loss Curve",
    )

    ax.legend(
        loc="lower left",
        ncol=3,
        labelspacing=0,
        columnspacing=0.5,
    )

    # ------------------  能量对比图 + 误差直方图  ------------------
    plt.subplot(2, 2, 2)
    ax = plt.gca()

    parity_plot(
        ax,
        energy_train[:, 1],
        energy_train[:, 0],
        "#2472A3",
        r"RMSE$_{\mathrm{train}}$",
    )

    if energy_test is not None:
        parity_plot(
            ax,
            energy_test[:, 1],
            energy_test[:, 0],
            "#EECA40",
            r"RMSE$_{\mathrm{test}}$",
        )

    xmin, xmax = ax.get_xlim()
    ax.plot([xmin, xmax], [xmin, xmax], c="black", lw=1)

    ax.set(
        xlim=[xmin, xmax],
        ylim=[xmin, xmax],
        xlabel="DFT energy (eV/atom)",
        ylabel="NEP energy (eV/atom)",
        title="(b) Energy",
    )

    ax.legend(loc="upper left")

    if hist:
        plot_hist(
            ax,
            energy_train[:, 1],
            energy_train[:, 0],
            r"$|E_{DFT}-E_{NEP}|$",
        )

    # ------------------  力对比图 + 误差直方图  ------------------
    plt.subplot(2, 2, 3)
    ax = plt.gca()

    parity_plot(
        ax,
        force_train[:, 3:].flatten(),
        force_train[:, :3].flatten(),
        "#2472A3",
        r"RMSE$_{\mathrm{train}}$",
    )

    if force_test is not None:
        parity_plot(
            ax,
            force_test[:, 3:].flatten(),
            force_test[:, :3].flatten(),
            "#EECA40",
            r"RMSE$_{\mathrm{test}}$",
        )

    xmin, xmax = ax.get_xlim()
    ax.plot([xmin, xmax], [xmin, xmax], c="black", lw=1)

    ax.set(
        xlim=[xmin, xmax],
        ylim=[xmin, xmax],
        xlabel=r"DFT force (eV/$\rm{\AA}$)",
        ylabel=r"NEP force (eV/$\rm{\AA}$)",
        title="(c) Force",
    )

    ax.legend(loc="upper left")

    if hist:
        plot_hist(
            ax,
            force_train[:, 3:].flatten(),
            force_train[:, :3].flatten(),
            r"$|F_{DFT}-F_{NEP}|$",
        )

    # ------------------  位力/应力对比图 + 误差直方图  ------------------
    plt.subplot(2, 2, 4)
    ax = plt.gca()

    if stress:
        virial_train = stress_train
        virial_test = stress_test
        axes_label = "stress (GPa)"
        title_label = "Stress"
    else:
        axes_label = "virial (eV/atom)"
        title_label = "Virial"

    parity_plot(
        ax,
        virial_train[:, 6:].flatten(),
        virial_train[:, :6].flatten(),
        "#2472A3",
        r"RMSE$_{\mathrm{train}}$",
    )

    if virial_test is not None:
        parity_plot(
            ax,
            virial_test[:, 6:].flatten(),
            virial_test[:, :6].flatten(),
            "#EECA40",
            r"RMSE$_{\mathrm{test}}$",
        )

    xmin, xmax = ax.get_xlim()
    ax.plot([xmin, xmax], [xmin, xmax], c="black", lw=1)

    ax.set(
        xlim=[xmin, xmax],
        ylim=[xmin, xmax],
        xlabel=f"DFT {axes_label}",
        ylabel=f"NEP {axes_label}",
        title=f"(d) {title_label}",
    )

    ax.legend(loc="upper left")

    if hist:
        plot_hist(
            ax,
            virial_train[:, 6:].flatten(),
            virial_train[:, :6].flatten(),
            r"$|V_{DFT}-V_{NEP}|$",
        )

    # ------------------  保存图片  ------------------
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.savefig("nep_loss_and_rmse.png", bbox_inches="tight")

    print("\nNEP train results figure is generated!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot NEP training results(loss evolution; parity plot of energy, force, virial, stress).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--stress",
        action="store_true",
        help="parity plot of stress or virial",
    )

    parser.add_argument(
        "--hist",
        action="store_true",
        help="Whether to plot the histogram of the error in subaxes",
    )

    args = parser.parse_args()

    plot_training_results(args.stress, args.hist)
