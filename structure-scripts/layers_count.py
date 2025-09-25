#!/usr/bin/env python3

"""统计原子层数目及每层原子数"""

import argparse

import numpy as np
import pandas as pd
from ase.io import read


def layers_count(
    structure_fn: str,
    precision: float = 0.001,
    axis: str = "z",
):
    """统计原子层数目及每层原子数"""

    if structure_fn.split(".")[-1] in ["lmp", "data", "lammps-data"]:
        atoms = read(structure_fn, format="lammps-data", sort_by_id=False)
    elif structure_fn.split(".")[-1] in ["xyz", "extxyz"]:
        atoms = read(structure_fn, format="extxyz")
    else:
        atoms = read(structure_fn)

    positions = atoms.get_positions(wrap=True)

    if axis == "z":
        positions_z = positions[:, 2]
    elif axis == "y":
        positions_z = positions[:, 1]
    elif axis == "x":
        positions_z = positions[:, 0]
    else:
        raise ValueError(f"Invalid axis: {axis}")

    positions_z_rounded = precision * np.round(positions_z / precision)

    positions_z_unique, natoms_layer_list = np.unique(
        positions_z_rounded, return_counts=True
    )

    count_layer = len(positions_z_unique)
    print(f"Atoms count: {len(atoms)}; Layers count: {count_layer}.\n")

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
        "structure_fn", nargs="?", default="POSCAR", help="structure filename"
    )

    parser.add_argument(
        "--axis",
        choices=["x", "y", "z"],
        default="z",
        metavar="STR",
        help="axis perpendicular to atomic layers (x, y, z)",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn
    axis = args.axis

    layers_count(structure_fn, axis=axis)
