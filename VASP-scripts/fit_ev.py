#!/usr/bin/env python3

"""
Birch-Murnaghan EOS 拟合

reference:
1. pymatgen.analysis.eos 模块
2. https://github.com/wangchr1617/learning/blob/main/scripts/birch_murnaghan_fitting.py
"""

import argparse

import numpy as np
import pandas as pd
from scipy.optimize import leastsq


def Birch_Murnaghan(params, volume):
    """Birch-Murnaghan 状态方程"""

    E0, V0, B0, B1 = params
    eta = (V0 / volume) ** (2 / 3)

    return E0 + (9 * V0 * B0 / 16) * (
        ((eta - 1) ** 3) * B1 + ((eta - 1) ** 2) * (6 - 4 * eta)
    )


def error(params, volume, energy):
    """误差函数"""

    return Birch_Murnaghan(params, volume) - energy


def initial_guess(volume, energy):
    """初始值猜测"""
    a, b, c = np.polyfit(volume, energy, 2)

    V0 = -b / (2 * a)
    E0 = a * V0**2 + b * V0 + c
    B0 = 2 * a * V0
    B1 = 4

    params_init = (E0, V0, B0, B1)

    return params_init


def fit(volume, energy):
    """EOS 拟合"""

    x0 = initial_guess(volume, energy)

    fitted_params = leastsq(
        func=error,
        x0=x0,
        args=(volume, energy),
    )

    E0, V0, B0, _ = fitted_params[0]
    coeffs_GPa = 1.60217662e-19 / (1e-30 * 1e9)
    B0 = round(B0 * coeffs_GPa, 1)

    print("EOS fitting results:\n")

    print(f"E0: {round(E0, 5): >10} eV/atom")
    print(f"V0: {round(V0, 5): >10} Å^3/atom")
    print(f"B0: {B0: >10} GPa")

    print("\nAssuming cubic box:")
    for i in range(1, 5):
        lattice_a = np.power(V0 * i, 1 / 3)
        print(f"If {i} atoms per cell, a = {lattice_a:.5f} Å")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Birch-Murnaghan EOS fitting",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "data_fn",
        type=str,
        nargs="?",
        default="ev.dat",
        help="File containing volume (Col2) and energy (Col3) data",
    )

    args = parser.parse_args()

    csv_fn = args.data_fn
    df = pd.read_csv(csv_fn, sep=None, engine="python")

    df_volume = df.iloc[:, 1]
    df_energy = df.iloc[:, 2]

    fit(df_volume, df_energy)
