#!/usr/bin/env python3

"""
通过 band.yaml 和 band.dat 文件绘制声子谱
band.yaml 文件由 phonopy -p band.conf -s 生成
band.dat 文件由 phonopy-bandplot --gnuplot > band.dat 生成
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import yaml

from spt.plot_params import set_plot_params


def plot_phonon_bandstructure(
    band_yaml_fn: str = "band.yaml",
    band_dat_fn: str = "band.dat",
    fig_fn: str = "band.png",
):

    with open(band_yaml_fn, "r") as f:
        yaml_data = yaml.safe_load(f)

    # K-PATH 路径对应的 x 轴坐标值；distances[0] 数目为 nqpoint
    distances = []
    distances.append([i["distance"] for i in yaml_data["phonon"]])

    # K-PATH labels
    labels = []
    for i, j in yaml_data["labels"]:
        labels.append(i)
        if len(labels) == len(yaml_data["labels"]):
            labels.append(j)
    labels = [r"$\Gamma$" if label == "G" else label for label in labels]

    # K-PATH label 对应的 x 轴坐标值；reference: phonopy/scripts/phonopy_bandplot.py 第 785 行附近
    end_points = [0]
    for nq in yaml_data["segment_nqpoint"]:
        end_points.append(nq + end_points[-1])
    # end_points 中的最后一个元素需要减 1
    end_points[-1] -= 1
    label_positions = [distances[0][i] for i in end_points]

    set_plot_params(roman_params=True, sci_params=True)
    fig, ax = plt.subplots(figsize=(10, 8))

    band_data = pd.read_csv(band_dat_fn, sep=" ", header=None, skiprows=2)

    # 以下方式可避免出现连接原点的线
    segment_nqpoint = yaml_data["segment_nqpoint"][0]
    for i in range(band_data.shape[0] // segment_nqpoint):
        ax.plot(
            band_data.iloc[i * segment_nqpoint : (i + 1) * segment_nqpoint, 0],
            band_data.iloc[i * segment_nqpoint : (i + 1) * segment_nqpoint, 1],
            color="#1f77b4",
        )

    ax.set_xticks(label_positions, labels, minor=False)
    ax.set_xlim(0, label_positions[-1])

    ax.set_ylabel("Frequency (THz)")

    ax.axhline(linewidth=2, color="black", linestyle="--")
    for label_positions in label_positions:
        ax.axvline(x=label_positions, linewidth=2, color="black")

    fig.savefig(fig_fn)

    print(f"\nPhonon band structure plot {fig_fn} is generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot phonon band structure with band.yaml and band.dat files."
    )

    parser.add_argument("band_yaml_fn", default="band.yaml", help="band.yaml filename")
    parser.add_argument("band_dat_fn", default="band.dat", help="band.dat filename")
    parser.add_argument(
        "fig_fn", default="band.png", help="phonon band structure figure filename"
    )

    args = parser.parse_args()

    plot_phonon_bandstructure(args.band_yaml_fn, args.band_dat_fn, args.fig_fn)
