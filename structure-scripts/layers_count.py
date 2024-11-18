#!/usr/bin/env python

"""统计原子层数目"""

import sys

import numpy as np
from ase.io import read


def layers_count(
    structure_fn: str,
    precision: float = 0.001,
):
    """统计原子层数目"""

    atoms = read(structure_fn)

    positions = atoms.positions
    positions_sorted = positions[positions[:, 2].argsort()]

    z_coords = positions_sorted[:, 2]
    z_coords_rounded = precision * np.round(z_coords / precision)

    z_unique = np.unique(z_coords_rounded)

    print(f"Number of atom layer : {len(z_unique)}")


if __name__ == "__main__":
    structure_fn = sys.argv[1]

    layers_count(structure_fn)
