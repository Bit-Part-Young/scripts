#!/usr/bin/env python3

"""
由 phonopy 命令行输出的 band.yaml 绘制声子谱，可叠加多个文件对比

band.yaml 由 `phonopy --config band.conf -s` 生成。
k 路径取自 band.yaml 本身；需要预设 k 路径、实验散点或多子图时用 plot_phonon_bs.py。
画图风格与 plot_phonon_bs.py 一致。

用法：
    plot_band_yaml.py band.yaml
    plot_band_yaml.py dft/band.yaml nep/band.yaml -l DFT NEP -o Nb_phonon.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

from plot_phonon_bs import DEFAULT_COLORS, DEFAULT_RC, tex_label

Loader = getattr(yaml, "CLoader", yaml.SafeLoader)


def read_band_yaml(band_yaml: str | Path) -> dict:
    """返回按段切分的 distances / frequencies 及高对称点标签"""

    with open(band_yaml) as f:
        data = yaml.load(f, Loader=Loader)

    distances = np.array([q["distance"] for q in data["phonon"]])
    frequencies = np.array([[b["frequency"] for b in q["band"]] for q in data["phonon"]])
    splits = np.cumsum(data["segment_nqpoint"])[:-1]

    return {
        "distances": np.split(distances, splits),
        "frequencies": np.split(frequencies, splits),
        "labels": data.get("labels"),  # [[start, end], ...]；未设 BAND_LABELS 时没有
    }


def get_xticks(band: dict) -> tuple[list[float], list[str]]:
    """每段起点 + 终点；不连续处合并为 'X|Y'"""

    positions = [d[0] for d in band["distances"]] + [band["distances"][-1][-1]]
    pairs = band["labels"]
    if pairs is None:
        return positions, [""] * len(positions)

    labels = [pairs[0][0]]
    for (_, end), (start, _) in zip(pairs[:-1], pairs[1:]):
        labels.append(end if end == start else f"{end}|{start}")
    labels.append(pairs[-1][1])

    return positions, ["|".join(tex_label(s) for s in label.split("|")) for label in labels]


def plot_band_yaml(
    band_yamls: list[str | Path],
    labels: list[str],
    output: str | Path,
    ylim: tuple[float, float] | None = None,
    figsize: tuple[float, float] = (3.4, 2.8),
):
    fig, ax = plt.subplots(figsize=figsize, layout="constrained")

    ref = None
    for i, (fn, label) in enumerate(zip(band_yamls, labels)):
        band = read_band_yaml(fn)
        # 各数据晶格常数略有差异，统一使用第一组的 x 坐标，保证高对称点对齐
        ref = ref or band
        nq = [len(d) for d in band["distances"]]
        if nq != [len(d) for d in ref["distances"]] or band["labels"] != ref["labels"]:
            raise ValueError(f"{fn}: k-path segments differ from {band_yamls[0]}")
        color = DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
        for j, (x, freq) in enumerate(zip(ref["distances"], band["frequencies"])):
            lines = ax.plot(x, freq, color=color, zorder=2 + len(band_yamls) - i)
            if j == 0:
                lines[0].set_label(label)

    positions, xticklabels = get_xticks(ref)
    for x in positions[1:-1]:
        ax.axvline(x, color="black", lw=0.5)
    ax.axhline(0, color="black", lw=0.5, ls="--")
    ax.set(
        xlim=(positions[0], positions[-1]),
        xticks=positions,
        xticklabels=xticklabels,
        ylabel="Frequency (THz)",
    )
    if ylim is not None:
        ax.set_ylim(ylim)
    if len(band_yamls) > 1:
        ax.legend()

    fig.savefig(output)
    print(f"Phonon band structure figure saved to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot phonon band structure from phonopy band.yaml file(s).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("band_yaml", nargs="+", help="band.yaml file(s)")
    parser.add_argument("-l", "--labels", nargs="+", help="legend labels (default: parent dir names)")
    parser.add_argument("-o", "--output_fn", default="band.png", help="output figure")
    parser.add_argument("--ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"))
    parser.add_argument("-s", "--style", nargs="+", help="matplotlib style file(s) or name(s)")
    args = parser.parse_args()

    labels = args.labels or [Path(fn).resolve().parent.name for fn in args.band_yaml]
    if len(labels) != len(args.band_yaml):
        parser.error("number of --labels must match band_yaml files")

    with plt.rc_context(DEFAULT_RC), plt.style.context(args.style or []):
        plot_band_yaml(args.band_yaml, labels, args.output_fn, ylim=args.ylim)
