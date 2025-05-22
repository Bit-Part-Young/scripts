"""
声子谱绘图（可绘制多组、可手动定义 K-path）

reference: https://github.com/wangchr1617/NEP_GT/blob/main/NEP_Phon/Fm-3m/phon.py
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seekpath
from ase.atoms import Atoms
from ase.io import read
from matplotlib.axes import Axes
from phonopy import load
from phonopy.api_phonopy import Phonopy
from phonopy.units import VaspToTHz

from spt.plot_params import set_plot_params


def get_kpoints(atoms: Atoms):
    """生成 K 点路径"""

    structure_tuple = (
        atoms.cell,
        atoms.get_scaled_positions(),
        atoms.numbers,
    )
    path = seekpath.get_explicit_k_path(structure_tuple)

    kpoints_rel, kpoints_lincoord, labels = (
        path["explicit_kpoints_rel"],
        path["explicit_kpoints_linearcoord"],
        path["explicit_kpoints_labels"],
    )
    labels = ["$\Gamma$" if label == "GAMMA" else label for label in labels]
    labels = [
        label.replace("_", "$_") + "$" if "_" in label else label for label in labels
    ]

    return kpoints_rel, kpoints_lincoord, labels


def get_manual_kpoints(high_symmetry_points, num=30):
    """
    基于手动定义的高对称点生成 K 点路径
    high_symmetry_points: 包含点坐标和标签的列表
    示例: BCC K-path
    high_symmetry_points = [
        {"label": "$\Gamma$", "coords": [0.0, 0.0, 0.0]},
        {"label": "H", "coords": [0.5, -0.5, 0.5]},
        {"label": "P", "coords": [0.25, 0.25, 0.25]},
        {"label": "$\Gamma$", "coords": [0.0, 0.0, 0.0]},
        {"label": "N", "coords": [0.0, 0.0, 0.5]},
    ]
    """

    kpoints_rel, kpoints_lincoord, labels = [], [], []
    current_distance = 0

    for i, point in enumerate(high_symmetry_points[:-1]):
        start, end = np.array(point["coords"]), np.array(
            high_symmetry_points[i + 1]["coords"]
        )
        kpoint_segment = np.linspace(start, end, num=num)
        kpoints_rel.extend(kpoint_segment)

        segment_length = np.linalg.norm(end - start)
        kpoints_lincoord.extend(
            np.linspace(current_distance, current_distance + segment_length, num=num)
        )

        current_distance += segment_length
        labels.extend([point["label"]] + [""] * (num - 1))
        if i == len(high_symmetry_points) - 2:
            labels[-1] = high_symmetry_points[i + 1]["label"]
        labels = ["$\Gamma$" if label == "GAMMA" else label for label in labels]

    return kpoints_rel, kpoints_lincoord, labels


def create_phonon_dataframe(phonon: Phonopy, kpoints_rel, kpoints_lincoord):
    """生成声子数据的DataFrame"""

    phonon.run_band_structure(
        paths=[kpoints_rel],
        with_eigenvectors=True,
        with_group_velocities=True,
    )

    bandstructure_dict = phonon.get_band_structure_dict()

    df = pd.DataFrame(bandstructure_dict["frequencies"][0], index=kpoints_lincoord)

    return df


def _plot_a_band(ax: Axes, df: pd.DataFrame, color: str, label: str, **kwargs):
    """绘制一组声子色散曲线"""

    for i, col in enumerate(df.columns):
        ax.plot(
            df.index,
            df[col],
            lw=3.0,
            color=color,
            label=label if i == 0 else None,
            **kwargs,
        )


def plot_phonon_dispersion(
    phonons: list[Phonopy],
    kpoints_rel,
    kpoints_lincoord,
    labels,
    fig_fn: str,
    **kwargs,
):
    """绘制（多组）声子色散曲线"""

    set_plot_params()

    # colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    colors = [
        "#bebebe",  # 灰色（用于不是很重要的数据）
        "#2b72bc",  # 蓝色（对比数据）
        "#fc0001",  # 红色（突出数据）
        # "#1f77b4",  # 蓝色 matplotlib
        # "#ff7f0e",  # 橙色 matplotlib
        # "#d62728",  # 红色 matplotlib
        # "#2ca02c",  # 绿色 matplotlib
    ]

    fig, ax = plt.subplots(figsize=(10, 8))

    for i, (label, phonon) in enumerate(phonons.items()):
        df = create_phonon_dataframe(phonon, kpoints_rel, kpoints_lincoord)
        _plot_a_band(ax, df, colors[i], label, **kwargs)

    df_path = pd.DataFrame(dict(labels=labels, positions=kpoints_lincoord))
    df_path = df_path[df_path.labels != ""]

    for xp in df_path.positions:
        ax.axvline(xp, c="black", lw=0.5)
    ax.axhline(y=0.0, c="black", lw=0.5)

    ax.set(
        xlabel="KPATH",
        ylabel="Frequency (THz)",
        xlim=(df.index.min(), df.index.max()),
        xticks=df_path.positions,
        xticklabels=df_path.labels,
    )
    ax.legend(loc="upper right")

    plt.tight_layout()

    fig.savefig(fig_fn)


if __name__ == "__main__":

    high_symmetry_points = [
        {"label": "$\Gamma$", "coords": [0.0, 0.0, 0.0]},
        {"label": "H", "coords": [0.5, -0.5, 0.5]},
        {"label": "P", "coords": [0.25, 0.25, 0.25]},
        {"label": "$\Gamma$", "coords": [0.0, 0.0, 0.0]},
        {"label": "N", "coords": [0.0, 0.0, 0.5]},
    ]

    kpoints_rel, kpoints_lincoord, labels = get_manual_kpoints(high_symmetry_points)

    # print(kpoints_rel)
    # print(kpoints_lincoord)
    # print(labels)

    phonon = load(
        "phonopy_disp.yaml",
        factor=VaspToTHz,
        is_nac=False,
        symprec=1e-4,
        force_sets_filename="FORCE_SETS",
    )

    phonons = {
        "FD": phonon,
    }

    plot_phonon_dispersion(
        phonons,
        kpoints_rel,
        kpoints_lincoord,
        labels,
        "phonon_disp.png",
    )
