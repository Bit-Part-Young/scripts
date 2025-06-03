#!/usr/bin/env python3

"""计算势函数预测与 DFT 计算的能量、力、应力/位力指标的 RMSE"""

import os
from glob import glob

import numpy as np


def calculate_rmse(data_fn):
    """计算 RMSE"""

    data = np.loadtxt(data_fn, ndmin=2)
    if "energy" in data_fn:
        pred, dft = data[:, 0], data[:, 1]
    elif "force" in data_fn:
        pred, dft = data[:, :3], data[:, 3:]
    elif "stress" in data_fn or "virial" in data_fn:
        pred, dft = data[:, :6], data[:, 6:]

    return np.sqrt(np.mean((pred - dft) ** 2))


def main():

    label_list = ["Energy", "Force", "Virial", "Stress"]
    unit_list = ["eV/atom", "eV/Å", "eV/atom", "GPa"]

    train_data_fn_list = [f"{label.lower()}_train.out" for label in label_list]
    test_data_fn_list = [f"{label.lower()}_test.out" for label in label_list]

    data_fn_list_list = [train_data_fn_list, test_data_fn_list]

    print("")
    for data_fn_list in data_fn_list_list:
        if os.path.exists(data_fn_list[0]):
            for data_fn, label, unit in zip(data_fn_list, label_list, unit_list):
                if os.path.exists(data_fn):
                    rmse = calculate_rmse(data_fn)
                    print(f"{label} train RMSE: {rmse:.6f} {unit}.")
            print("")


if __name__ == "__main__":
    main()
