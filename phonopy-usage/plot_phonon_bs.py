#!/usr/bin/env python3

"""
声子谱绘制

reference: https://github.com/JaGeo/mace-mp-03b-phonon-benchmark/blob/main/functions/phonon.py
"""

import argparse
import copy

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from phonopy import load
from phonopy.api_phonopy import Phonopy
from phonopy.phonon.band_structure import get_band_qpoints, get_band_qpoints_by_seekpath
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.units import VaspToTHz
from pymatgen.core.structure import Structure
from pymatgen.io.phonopy import get_pmg_structure
from pymatgen.symmetry.kpath import KPathSeek
from spt.plot_params import set_plot_params


def load_phonopy(
    phonopy_yaml_fn: str,
    symprec: float = 1e-4,
    is_nac: bool = False,
    force_sets_fn: str | None = None,
    force_constants_fn: str | None = None,
) -> Phonopy:

    if force_sets_fn is not None:
        phonon = load(
            phonopy_yaml=phonopy_yaml_fn,
            factor=VaspToTHz,
            is_nac=is_nac,
            symprec=symprec,
            force_sets_filename=force_sets_fn,
        )
    elif force_constants_fn is not None:
        phonon = load(
            phonopy_yaml=phonopy_yaml_fn,
            factor=VaspToTHz,
            is_nac=is_nac,
            symprec=symprec,
            force_constants_filename=force_constants_fn,
        )
    else:
        phonon = load(
            phonopy_yaml=phonopy_yaml_fn,
            factor=VaspToTHz,
            is_nac=is_nac,
            symprec=symprec,
            force_sets_filename=phonopy_yaml_fn,
        )

    return phonon


# 常见晶体结构的标准 k-path（原胞倒格矢分数坐标，Setyawan-Curtarolo 约定）
KPATH_PRESETS: dict[str, dict] = {
    "bcc": {
        "kpoints": {
            "$\\Gamma$": [0, 0, 0],
            "H": [1 / 2, -1 / 2, 1 / 2],
            "P": [1 / 4, 1 / 4, 1 / 4],
            "N": [0, 0, 1 / 2],
        },
        "path": [
            ["$\\Gamma$", "H", "N", "$\\Gamma$", "P"],
            ["P", "H"],
            ["P", "N"],
        ],
    },
    "bcc_partial": {
        "kpoints": {
            "$\\Gamma$": [0, 0, 0],
            "H": [1 / 2, -1 / 2, 1 / 2],
            "P": [1 / 4, 1 / 4, 1 / 4],
            "N": [0, 0, 1 / 2],
        },
        "path": [
            ["$\\Gamma$", "H", "P", "$\\Gamma$", "N"],
        ],
    },
    "fcc": {
        "kpoints": {
            "$\\Gamma$": [0, 0, 0],
            "X": [1 / 2, 0, 1 / 2],
            "W": [1 / 2, 1 / 4, 3 / 4],
            "K": [3 / 8, 3 / 8, 3 / 4],
            "L": [1 / 2, 1 / 2, 1 / 2],
            "U": [5 / 8, 1 / 4, 5 / 8],
        },
        "path": [
            ["$\\Gamma$", "X", "W", "K", "$\\Gamma$", "L"],
            ["L", "U", "W", "L", "K"],
            ["U", "X"],
        ],
    },
    "fcc_partial": {
        "kpoints": {
            "$\\Gamma$": [0, 0, 0],
            "X": [1 / 2, 0, 1 / 2],
            "U": [5 / 8, 1 / 4, 5 / 8],
            "K": [3 / 8, 3 / 8, 3 / 4],
            "L": [1 / 2, 1 / 2, 1 / 2],
        },
        "path": [
            ["$\\Gamma$", "X", "U"],
            ["K", "$\\Gamma$", "L"],
        ],
        "label_map": {"U | K": "K"},  # 用于将 FCC 的 "U | K" xtick label 替换为 "K"
    },
    "hcp": {
        "kpoints": {
            "$\\Gamma$": [0, 0, 0],
            "A": [0, 0, 1 / 2],
            "M": [1 / 2, 0, 0],
            "L": [1 / 2, 0, 1 / 2],
            "K": [1 / 3, 1 / 3, 0],
            "H": [1 / 3, 1 / 3, 1 / 2],
        },
        "path": [
            ["$\\Gamma$", "M", "K", "$\\Gamma$", "A", "L", "H", "A"],
            ["L", "M"],
            ["K", "H"],
        ],
    },
    "hcp_partial": {
        "kpoints": {
            "$\\Gamma$": [0, 0, 0],
            "A": [0, 0, 1 / 2],
            "M": [1 / 2, 0, 0],
            "K": [1 / 3, 1 / 3, 0],
        },
        "path": [
            ["$\\Gamma$", "K", "M", "$\\Gamma$", "A"],
        ],
    },
}


def get_kpath_preset(structure_type: str) -> tuple[list, list[str], list[str]]:
    """根据晶体结构类型返回预定义的 k-path、标签和 xtick 标签。

    phonopy 对不连续路径段共享同一 distance 端点（无间隔），因此需要在断点处
    将前后段的标签合并为 "X | Y" 形式，保证标签数与 xtick 位置数一致。

    Returns:
        path: 分段 k-path 坐标列表，供 get_band_qpoints 使用
        labels: 全部高对称点标签的平展列表，供 phonon.run_band_structure 使用
        xtick_labels: 用于 ax.set_xticklabels 的最终标签列表
    """

    key = structure_type.lower()
    if key not in KPATH_PRESETS:
        raise ValueError(
            f"Unknown structure type: '{structure_type}'. "
            f"Available: {list(KPATH_PRESETS.keys())}"
        )

    preset = KPATH_PRESETS[key]
    kpoints = preset["kpoints"]
    label_path = preset["path"]

    path = []
    labels: list[str] = []
    xtick_labels: list[str] = list(label_path[0])

    for seg_idx, segment in enumerate(label_path):
        coords = [kpoints[label] for label in segment]
        path.append(coords)
        labels.extend(segment)

        if seg_idx > 0:
            prev_end = xtick_labels[-1]
            curr_start = segment[0]
            if prev_end != curr_start:
                xtick_labels[-1] = f"{prev_end} | {curr_start}"
            xtick_labels.extend(segment[1:])

    # 用于将 FCC 的 "U | K" xtick label 替换为 "K"
    label_map = preset.get("label_map", {})
    xtick_labels = [label_map.get(lbl, lbl) for lbl in xtick_labels]

    return path, labels, xtick_labels


def get_kpath(primitive: PhonopyAtoms) -> list:
    structure: Structure = get_pmg_structure(primitive)

    kpath_high_symmetry = KPathSeek(structure=structure, symprec=1e-4)
    kpath = kpath_high_symmetry.kpath

    path = copy.deepcopy(kpath["path"])
    for idx, labelset in enumerate(kpath["path"]):
        for i, label in enumerate(labelset):
            path[idx][i] = kpath["kpoints"][label]

    return path


def get_label_and_connection(primitive: PhonopyAtoms, npoints: int = 101) -> tuple:
    bands, labels_for_plot, connections = get_band_qpoints_by_seekpath(
        primitive=primitive,
        npoints=npoints,
        is_const_interval=True,
    )

    return labels_for_plot, connections


def run_bands_structure_dict(
    phonon: Phonopy, path, labels: list, npoints: int = 101
) -> dict:
    band_qpoints = get_band_qpoints(band_paths=path, npoints=npoints)

    phonon.run_band_structure(
        paths=band_qpoints,
        labels=labels,
        is_band_connection=True,
    )

    bandstructure_dict = phonon.get_band_structure_dict()

    return bandstructure_dict


def add_vertical_lines_and_commensurate_points(
    ax: Axes, distances: np.ndarray, comm_points
):
    """add vertical lines and commensurate points to the plot"""

    xticks = []
    num_paths = distances.shape[0]
    for path in range(num_paths):
        start, end = distances[path][0], distances[path][-1]
        ax.axvline(x=start, color="black", linestyle="--", linewidth=0.5)
        ax.axvline(x=end, color="black", linestyle="--", linewidth=0.5)
        xticks.extend([start, end])
        if comm_points is not None and len(comm_points) > 0:
            for point in comm_points[path]:
                ax.plot(distances[path][point], 0, color="red", marker="D")

    return sorted(set(xticks))


def create_xtick_labels(x_labels: list[str], connections: list[bool]) -> list[str]:
    """create xtick labels for the phonon band structure plot
    Depending on wehter the x-path is connected or not, the labels have to repeat or change
    An example:
    x_labels = [Gamma,X,K,L, Gamma]
    connection = [True, True, False]
    '--------'---------'-----------'
    G        X        K| L         G

    x_labels: list of labels for example: [Gamma,X,K, Gamma]
    connection: list with True and False label if x-path is connected or not. for example: [True, True, False, True]
    """

    xtick_labels = []
    if False not in connections:
        return x_labels
    else:
        xtick_labels.append(x_labels[0])
        count = 1

        for connection in connections:
            if count >= len(x_labels) - 1:
                xtick_labels.append(x_labels[-1])
                break
            if connection:
                xtick_labels.append(x_labels[count])
                count += 1
            else:
                xtick_labels.append(str(x_labels[count]) + " | " + x_labels[count + 1])
                count += 2

    return xtick_labels


def plot_phonon_bandstructure(
    phonon_distance, phonon_frequencies, xtick_labels, figure_fn: str = "phonon_bs.png"
):

    set_plot_params(roman_params=True, sci_params=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    num_paths, num_kpoints, num_bands = np.array(phonon_frequencies).shape

    for path in range(num_paths):
        for band in range(num_bands):
            ax.plot(
                phonon_distance[path],
                np.array(phonon_frequencies)[path, :, band],
                color="black",
            )

    # 添加垂直线及其对应的 x 轴刻度、标签
    commensurate_points = None
    default_comm_points = (
        commensurate_points
        if commensurate_points is not None
        else [[]] * np.array(phonon_distance).shape[0]
    )

    xticks = add_vertical_lines_and_commensurate_points(
        ax, np.array(phonon_distance), default_comm_points
    )

    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)

    ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)
    ax.set_ylabel("Frequency [THz]")
    # ax.set_xlabel("Wave vector")
    ax.set_xlim(0, np.array(phonon_distance)[-1, -1])
    # ax.legend(loc="lower right")

    fig.savefig(figure_fn)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Phonon bandstructure plot with phonopy.",
    )

    parser.add_argument(
        "-i",
        "--input_fn",
        nargs="+",
        metavar="FILE",
        help="1 or 2 phonon data files; eg. phonopy*.yaml, FORCE_SETS, FORCE_CONSTANTS",
    )
    parser.add_argument(
        "-o",
        "--output_fn",
        default="phonon_bs.png",
        metavar="FILE",
        help="output figure filename (default: phonon_bs.png)",
    )
    parser.add_argument(
        "-pk",
        "--preset_kpath",
        default=None,
        choices=list(KPATH_PRESETS.keys()),
        help="use preset/predefined k-path for the given bcc/fcc/hcp structure; "
        "eg. 'bcc' is the full k-path, 'bcc_partial' is the partial k-path in literature "
        "(default: auto-detect via SeeK-path)",
    )

    args = parser.parse_args()
    input_fn = args.input_fn
    figure_fn = args.output_fn

    phonon_yaml_fn = input_fn[0]
    if len(input_fn) == 1:
        phonon = load_phonopy(phonon_yaml_fn)
    elif len(input_fn) == 2:
        if input_fn[1] == "FORCE_SETS":
            phonon = load_phonopy(
                phonopy_yaml_fn=phonon_yaml_fn,
                force_sets_fn=input_fn[1],
            )
        elif input_fn[1] == "FORCE_CONSTANTS":
            phonon = load_phonopy(
                phonopy_yaml_fn=phonon_yaml_fn,
                force_constants_fn=input_fn[1],
            )

    if args.structure is not None:
        path, labels, xtick_labels = get_kpath_preset(args.structure)
    else:
        path = get_kpath(phonon.primitive)
        labels, connections = get_label_and_connection(phonon.primitive)
        xtick_labels = create_xtick_labels(labels, connections)

    bandstructure_dict = run_bands_structure_dict(
        phonon=phonon, path=path, labels=labels
    )
    phonon_distance = bandstructure_dict["distances"]
    phonon_frequencies = bandstructure_dict["frequencies"]

    plot_phonon_bandstructure(
        phonon_distance, phonon_frequencies, xtick_labels, figure_fn
    )
