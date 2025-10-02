#!/usr/bin/env python3

"""将 EOS 中的 e-v 数据转换为 deltaE-v 数据（用于绘制 deltaE-v EOS 曲线）"""

import argparse

import numpy as np
import pandas as pd
from scipy.optimize import leastsq


def Birch_Murnaghan(params, volumes):
    """Birch-Murnaghan EOS 状态方程"""

    E0, V0, B0, B1 = params
    eta = (V0 / volumes) ** (2 / 3)

    return E0 + (9 * V0 * B0 / 16) * (
        ((eta - 1) ** 3) * B1 + ((eta - 1) ** 2) * (6 - 4 * eta)
    )


def error(params, volumes, energies):
    """误差函数"""

    return Birch_Murnaghan(params, volumes) - energies


def initial_guess(volumes, energies):
    """初始值猜测"""
    a, b, c = np.polyfit(volumes, energies, 2)

    V0 = -b / (2 * a)
    E0 = a * V0**2 + b * V0 + c
    B0 = 2 * a * V0
    B1 = 4

    params_init = (E0, V0, B0, B1)

    return params_init


def fit(volumes, energies):
    """EOS 拟合"""

    x0 = initial_guess(volumes, energies)

    fitted_params = leastsq(func=error, x0=x0, args=(volumes, energies))

    E0, V0, B0, _ = fitted_params[0]
    coeffs_GPa = 1.602176634e-19 / (1e-30 * 1e9)
    B0 = round(B0 * coeffs_GPa, 1)

    return round(E0, 5)


def deltaE_generation(volumes: pd.Series, energies: pd.Series, output_data_fn: str):
    """生成 deltaE 数据"""

    E0 = fit(volumes, energies)

    df_deltaE = energies - E0

    df = pd.DataFrame({"volume": volumes, "deltaE": df_deltaE}).round(5)

    df.to_csv(output_data_fn, index=False, sep=" ")

    print(f"\n{output_data_fn} is generated.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Transfer e-v data in EOS to deltaE-v.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input_data_fn", default="ev.dat", help="input EOS data filename"
    )
    parser.add_argument(
        "cols",
        type=int,
        nargs=2,
        default=[2, 3],
        help="Column numbers for volume and energy data",
    )
    parser.add_argument(
        "output_data_fn", default="deltaE_v.dat", help="output deltaE data filename"
    )

    args = parser.parse_args()

    input_data_fn = args.input_data_fn
    cols = args.cols
    output_data_fn = args.output_data_fn

    df = pd.read_csv(input_data_fn, sep=None, engine="python")
    df_volume = df.iloc[:, cols[0] - 1]
    df_energy = df.iloc[:, cols[1] - 1]

    deltaE_generation(df_volume, df_energy, output_data_fn)
