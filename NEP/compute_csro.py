#!/usr/bin/env python3

"""
计算 BCC 结构 MC/MD 模拟轨迹的 1NN CSRO / WCP 参数

reference: https://github.com/killiansheriff/WarrenCowleyParameters
"""

import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

import argparse

import pandas as pd
import WarrenCowleyParameters as wc
from ovito.io import import_file


def compute_csro(trajectory_fn: str = "dump.xyz", save_output: bool | None = None):
    """计算 BCC 结构的 1NN CSRO / WCP 参数"""

    # trajectory (GPUMD / LAMMPS)
    pipeline = import_file(trajectory_fn)

    nframes = pipeline.num_frames
    print(f"Number of frames: {nframes}.\n")

    # BCC 1NN 2NN
    modifier = wc.WarrenCowleyParameters(nneigh=[0, 8, 14], only_selected=False)
    # FCC 1NN 2NN
    # modifier = wc.WarrenCowleyParameters(nneigh=[0, 12, 18], only_selected=False)

    pipeline.modifiers.append(modifier)

    # 目前是处理所有帧
    wc_1NN_list = []
    for i in range(nframes):

        data = pipeline.compute(i)

        wc_for_shells = data.attributes["Warren-Cowley parameters"]

        print(f"\nNo. {i} processed.\n")

        wc_list = data.attributes["Warren-Cowley parameters by particle name"]

        wc_1NN_dict = wc_list[0]
        wc_1NN_list.append(wc_1NN_dict)

    df = pd.DataFrame(wc_1NN_list).round(5)
    print(df)

    if save_output:
        df.to_csv("wc_1NN.csv", index=False, sep=" ")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute CSRO / Warren-Cowley parameters for given trajectory in MC/MD simulation.",
        epilog="Author: SLY.",
    )

    parser.add_argument("trajectory_fn", default="dump.xyz", help="trajectory filename")
    parser.add_argument("-o", action="store_true", help="save 1NN CSRO / WCP data file")

    args = parser.parse_args()

    compute_csro(args.trajectory_fn, args.o)
