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
from phonopy.units import VaspToTHz
from pymatgen.io.phonopy import get_pmg_structure
from pymatgen.symmetry.kpath import KPathSeek

from spt.plot_params import set_plot_params


def load_phonopy(
    phonopy_yaml_filename: str,
    symprec: float = 1e-4,
    is_nac: bool = False,
    force_sets_filename: str | None = None,
    force_constants_filename: str | None = None,
) -> Phonopy:

    if force_sets_filename is not None:
        phonon = load(
            phonopy_yaml=phonopy_yaml_filename,
            factor=VaspToTHz,
            is_nac=is_nac,
            symprec=symprec,
            force_sets_filename=force_sets_filename,
        )
    elif force_constants_filename is not None:
        phonon = load(
            phonopy_yaml=phonopy_yaml_filename,
            factor=VaspToTHz,
            is_nac=is_nac,
            symprec=symprec,
            force_constants_filename=force_constants_filename,
        )
    else:
        phonon = load(
            phonopy_yaml=phonopy_yaml_filename,
            factor=VaspToTHz,
            is_nac=is_nac,
            symprec=symprec,
            force_sets_filename=phonopy_yaml_filename,
        )

    return phonon


def get_kpath(primitive) -> list:
    structure_pmg = get_pmg_structure(primitive)

    highsymmkpath = KPathSeek(structure=structure_pmg, symprec=1e-4)
    kpath = highsymmkpath.kpath
    # print(kpath)

    path = copy.deepcopy(kpath["path"])
    for idx, labelset in enumerate(kpath["path"]):
        for i, label in enumerate(labelset):
            path[idx][i] = kpath["kpoints"][label]

    # print(path)

    return path


def get_label_and_connection(primitive, npoints: int = 101) -> tuple:
    bands, labels_for_plot, connections = get_band_qpoints_by_seekpath(
        primitive=primitive,
        npoints=npoints,
        is_const_interval=True,
    )

    return labels_for_plot, connections


def run_bands_structure_dict(phonon: Phonopy, path, labels, npoints: int = 101) -> dict:
    band_qpoints = get_band_qpoints(band_paths=path, npoints=npoints)

    phonon.run_band_structure(
        paths=band_qpoints,
        labels=labels,
        is_band_connection=True,
    )

    bandstructure_dict = phonon.get_band_structure_dict()
    # print(bandstructure_dict)

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


def plot_phonon_bs(phonon_distance, phonon_frequencies, xtick_labels):

    set_plot_params()

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

    # additional plot settings
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)
    ax.set_ylabel("Frequency [THz]")
    ax.set_xlabel("Wave vector")
    ax.set_xlim(0, np.array(phonon_distance)[-1, -1])
    # ax.legend(loc="lower right")

    fig.savefig("phonon_bands.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Phonon bandstructure plot with phonopy.",
    )

    parser.add_argument(
        "-i",
        "--input_fn",
        nargs="+",
        type=str,
        help="1 or 2 phonopy input files",
    )

    args = parser.parse_args()
    input_fn = args.input_fn

    phonon_yaml_filename = input_fn[0]
    if len(input_fn) == 1:
        phonon = load_phonopy(phonon_yaml_filename)
    elif len(input_fn) == 2:
        if input_fn[1] == "FORCE_SETS":
            phonon = load_phonopy(
                phonopy_yaml_filename=phonon_yaml_filename,
                force_sets_filename=input_fn[1],
            )
        elif input_fn[1] == "FORCE_CONSTANTS":
            phonon = load_phonopy(
                phonopy_yaml_filename=phonon_yaml_filename,
                force_constants_filename=input_fn[1],
            )

    path = get_kpath(phonon.primitive)
    labels_for_plot, connections = get_label_and_connection(phonon.primitive)

    bandstructure_dict = run_bands_structure_dict(
        phonon=phonon,
        path=path,
        labels=labels_for_plot,
    )
    phonon_distance = bandstructure_dict["distances"]
    phonon_frequencies = bandstructure_dict["frequencies"]

    xtick_labels = create_xtick_labels(labels_for_plot, connections)
    # print(xtick_labels)

    plot_phonon_bs(phonon_distance, phonon_frequencies, xtick_labels)
