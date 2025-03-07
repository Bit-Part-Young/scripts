#!/usr/bin/env python3

"""计算势函数预测与DFT 计算的能量、力、应力、位力指标的 RMSE"""

import os

import numpy as np


def calculate_rmse(data_fn):
    """计算 RMSE"""

    data = np.loadtxt(data_fn)
    if data_fn == "energy_train.out":
        pred, dft = data[:, 0], data[:, 1]
    elif data_fn == "force_train.out":
        pred, dft = data[:, :3], data[:, 3:]
    elif data_fn in ["virial_train.out", "stress_train.out"]:
        pred, dft = data[:, :6], data[:, 6:]

    return np.sqrt(np.mean((pred - dft) ** 2))


def main():

    data_fn_list = [
        "energy_train.out",
        "force_train.out",
        "virial_train.out",
        "stress_train.out",
    ]

    label_list = ["Energy", "Force", "Virial", "Stress"]
    unit_list = ["eV/atom", "eV/Å", "eV/atom", "GPa"]

    for data_fn, label, unit in zip(data_fn_list, label_list, unit_list):
        if os.path.exists(data_fn):
            rmse = calculate_rmse(data_fn)
            print(f"{label} RMSE: {rmse:.6f} {unit}.")


if __name__ == "__main__":
    main()
