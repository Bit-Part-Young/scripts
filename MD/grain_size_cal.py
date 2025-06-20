#!/usr/bin/env python3

"""计算由 atomsk 生成的多晶的平均晶粒尺寸"""

from glob import glob

import numpy as np


def grain_size_cal():
    """计算由 atomsk 生成的多晶的平均晶粒尺寸"""

    data_fn = glob("*_id-size.txt")[0]
    data = np.loadtxt(data_fn)

    grain_volume = data[:, 2]
    grain_size = np.power(grain_volume * 3 / (4 * np.pi), 1 / 3)
    grain_size_mean = grain_size.mean()
    grain_size_std = grain_size.std()

    print(f"Average Grain Size: {grain_size_mean*2:.2f} ± {grain_size_std*2:.2f} Å.")
    print(
        f"Average Grain Size: {grain_size_mean*2/10:.2f} ± {grain_size_std*2/10:.2f} nm."
    )


if __name__ == "__main__":
    grain_size_cal()
