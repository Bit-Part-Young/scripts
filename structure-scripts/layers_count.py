#!/usr/bin/env python

"""
统计原子层数目及每层原子数

Usage:
    layers_count.py POSCAR
"""

import sys

import numpy as np
import pandas as pd
from ase.io import read


def layers_count(
    structure_fn: str,
    precision: float = 0.001,
):
    """
    统计原子层数目及每层原子数
    """

    atoms = read(structure_fn)

    positions = atoms.positions

    z_coords = positions[:, 2]
    z_coords_rounded = precision * np.round(z_coords / precision)

    z_unique, layer_natoms_list = np.unique(
        z_coords_rounded,
        return_counts=True,
    )

    count_layer = len(z_unique)
    print(f"Layers conut: {count_layer}.\n")

    index = [f"Layer{i}" for i in range(1, count_layer + 1)]
    df = pd.DataFrame(
        {"natoms": layer_natoms_list},
        index=index,
    )

    print(df)


if __name__ == "__main__":
    structure_fn = sys.argv[1]

    layers_count(structure_fn)
