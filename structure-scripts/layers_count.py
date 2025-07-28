#!/usr/bin/env python3

"""统计原子层数目及每层原子数"""

import argparse

import numpy as np
import pandas as pd
from ase.io import read


def layers_count(
    structure_fn: str,
    precision: float = 0.001,
):
    """统计原子层数目及每层原子数"""

    if structure_fn.split(".")[-1] in ["lmp", "lammps-data"]:
        atoms = read(structure_fn, format="lammps-data")
    elif structure_fn.split(".")[-1] in ["xyz", "extxyz"]:
        atoms = read(structure_fn, format="extxyz")
    else:
        atoms = read(structure_fn)

    positions = atoms.positions

    positions_z = positions[:, 2]
    positions_z_rounded = precision * np.round(positions_z / precision)

    positions_z_unique, natoms_layer_list = np.unique(
        positions_z_rounded, return_counts=True
    )

    count_layer = len(positions_z_unique)
    print(f"Atoms count: {atoms.num_atoms}; Layers count: {count_layer}.\n")

    index = [f"layer{i}" for i in range(1, count_layer + 1)]
    df = pd.DataFrame({"natoms": natoms_layer_list}, index=index)

    print(df)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Count atomic layers and the number of atoms in each layer.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn",
        nargs="?",
        default="POSCAR",
        metavar="structure_fn",
        help="structure filename",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn

    layers_count(structure_fn)
