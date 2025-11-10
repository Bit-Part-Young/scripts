#!/usr/bin/env python3

"""
通过 band.yaml 文件绘制声子谱
reference: https://github.com/warda-rahim/phononplotter/blob/master/phonon.py
"""

import argparse
import yaml
import matplotlib.pyplot as plt
from spt.plot_params import set_plot_params


def plot_phonon_bandstructure(band_yaml_fn: str, fig_fn: str = "band.png"):

    with open(band_yaml_fn, "r") as f:
        data = yaml.safe_load(f)

    # K-PATH 路径对应的 x 轴坐标值；distances[0] 数目为 nqpoint
    distances = []
    distances.append([i["distance"] for i in data["phonon"]])

    # K-PATH 路径对应的频率值（本征值）；数目为 nqpoint * (natoms * 3)
    eigenvalues = []
    for point in data["phonon"]:
        eigenvalues.append([e["frequency"] for e in point["band"]])

    # K-PATH labels
    labels = []
    for i, j in data["labels"]:
        labels.append(i)
        if len(labels) == len(data["labels"]):
            labels.append(j)
    labels = [r"$\Gamma$" if label == "G" else label for label in labels]

    # K-PATH label 对应的 x 轴坐标值
    label_positions = []
    step = data["segment_nqpoint"][0]
    for i in range(0, len(distances[0]), step - 1):
        label_positions.append(distances[0][i])

    set_plot_params(roman_params=True, sci_params=True)
    fig, ax = plt.subplots(figsize=(10, 8))

    ax.plot(distances[0], eigenvalues, color="#1f77b4")

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
        description="Plot phonon band structure with band.yaml file."
    )

    parser.add_argument("band_yaml_fn", default="band.yaml", help="band.yaml filename")
    parser.add_argument(
        "fig_fn", default="band.png", help="phonon band structure figure filename"
    )

    args = parser.parse_args()

    plot_phonon_bandstructure(args.band_yaml_fn, args.fig_fn)
