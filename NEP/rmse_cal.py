#!/usr/bin/env python3

"""计算 NEP 势函数预测与 DFT 计算的能量、力、应力/位力指标的 RMSE MAE R2"""

import argparse
import os

import numpy as np


def statistics_calc(data_fn, statistic_type: str = "rmse"):
    """计算 RMSE MAE R2"""

    data = np.loadtxt(data_fn, ndmin=2)
    if "energy" in data_fn:
        pred, dft = data[:, 0], data[:, 1]
    elif "force" in data_fn:
        pred, dft = data[:, :3], data[:, 3:]
    elif "stress" in data_fn or "virial" in data_fn:
        pred, dft = data[:, :6], data[:, 6:]

    if statistic_type == "rmse":
        return np.sqrt(np.mean((pred - dft) ** 2))
    elif statistic_type == "mae":
        return np.mean(np.abs(pred - dft))
    elif statistic_type == "r2":
        return 1 - np.sum((pred - dft) ** 2) / np.sum((dft - np.mean(dft)) ** 2)


def statistics_info(statistic_type: str = "rmse"):
    label_list = ["Energy", "Force", "Virial", "Stress"]
    unit_list = ["eV/atom", "eV/Å", "eV/atom", "GPa"]

    train_data_fn_list = [f"{label.lower()}_train.out" for label in label_list]
    test_data_fn_list = [f"{label.lower()}_test.out" for label in label_list]

    data_fn_list_list = [train_data_fn_list, test_data_fn_list]

    print()
    for data_fn_list in data_fn_list_list:
        if os.path.exists(data_fn_list[0]):
            for data_fn, label, unit in zip(data_fn_list, label_list, unit_list):
                if os.path.exists(data_fn):
                    result = statistics_calc(data_fn, statistic_type)
                    tag = "train"
                    if "test" in data_fn:
                        tag = "test"
                    print(
                        f"{label} {tag} {statistic_type.upper()}: {result:.5f} {unit}."
                    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate RMSE MAE R2 of NEP prediction and DFT values for energy, force, stress/virial.",
        epilog="Author: SLY.",
    )

    args = parser.parse_args()

    statistics_info(statistic_type="rmse")
    statistics_info(statistic_type="mae")
    statistics_info(statistic_type="r2")
