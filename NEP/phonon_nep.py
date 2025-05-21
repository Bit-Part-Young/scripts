#!/usr/bin/env python3

"""
NEP 势函数声子谱计算（使用 calorine + phonopy）

Original Author: Ke Xu
"""

import argparse
import os

import numpy as np
from ase.atoms import Atoms
from ase.io import read
from calorine.calculators import CPUNEP
from calorine.tools import get_force_constants, relax_structure
from pymatgen.core.structure import Structure
from seekpath import get_explicit_k_path


def get_primitive(atoms: Atoms):
    """获取原胞"""

    structures = Structure.from_ase_atoms(atoms)

    primitive = structures.get_primitive_structure()

    return primitive.to_ase_atoms()


def get_kpath(atoms: Atoms, custom_path=None):

    if custom_path is not None:
        # 自定义 k 点路径插值计算
        explicit_kpoints_rel = []
        explicit_labels = []
        num_points_per_segment = 68

        for start_label, end_label in custom_path:
            # 线性插值生成 k 点路径
            start = np.array(point_coords[start_label])
            end = np.array(point_coords[end_label])
            kpoints = [
                start + (end - start) * i / (num_points_per_segment - 1)
                for i in range(num_points_per_segment)
            ]
            explicit_kpoints_rel.extend(kpoints)
            # 路径标签处理
            if not explicit_labels or explicit_labels[-1] != start_label:
                explicit_labels.append(start_label)
            explicit_labels.extend([""] * (num_points_per_segment - 2))
            explicit_labels.append(end_label)
        explicit_kpoints_rel = np.array(explicit_kpoints_rel)
        # Create linear coordinates for plotting
        explicit_kpoints_linearcoord = np.linspace(0, 1, len(explicit_kpoints_rel))
        path = {
            "explicit_kpoints_rel": explicit_kpoints_rel,
            "explicit_kpoints_labels": explicit_labels,
            "explicit_kpoints_linearcoord": explicit_kpoints_linearcoord,
        }
    else:
        structure_tuple = (
            atoms.cell,
            atoms.get_scaled_positions(),
            atoms.numbers,
        )
        path = get_explicit_k_path(structure_tuple, reference_distance=0.02)
        point_coords = path["point_coords"]

    return path


def get_phonon_dispersion(
    atoms: Atoms,
    model_fn,
    cellsize=4,
    custom_path=None,
    phonopy_yaml_fn="phonopy_params.yaml",
):
    """计算声子色散"""

    calculator = CPUNEP(model_fn)

    atoms = get_primitive(atoms)

    atoms.calc = calculator

    relax_structure(
        structure=atoms,
        fmax=0.00001,
        minimizer="bfgs",
        constant_volume=False,
        constant_cell=False,
    )

    path = get_kpath(atoms, custom_path=custom_path)

    # 声子谱计算
    phonon = get_force_constants(
        structure=atoms,
        calculator=calculator,
        supercell_matrix=[cellsize, cellsize, cellsize],
    )
    phonon.run_band_structure([path["explicit_kpoints_rel"]])

    # 将 phonopy 的相关参数写入文件
    phonon.save(filename=phonopy_yaml_fn)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="NEP phonon calculation with calorine & phonopy."
    )

    parser.add_argument("structure_fn", type=str, help="structure filename")

    args = parser.parse_args()

    atoms = read(args.structure_fn)
    cal_type = "NEP"
    element = "-".join(list(set(atoms.get_chemical_symbols())))

    folder = f"{element}-{cal_type}"
    os.makedirs(folder, exist_ok=True)

    get_phonon_dispersion(
        atoms,
        model_fn="nep.txt",
        cellsize=4,
        phonopy_yaml_fn=f"{folder}/phonopy_params.yaml",
    )
