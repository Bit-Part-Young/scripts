#!/usr/bin/env python

"""
统计原子层间距变化
仅限表面构型弛豫前后的原子层间距变化
"""

import sys

import numpy as np
from ase.io import read


def get_distance(
    structure_fn: str,
    precision: float = 0.001,
):
    """统计原子层间距"""

    atoms = read(structure_fn)

    positions = atoms.positions
    positions_sorted = positions[positions[:, 2].argsort()]

    z_coords = positions_sorted[:, 2]
    z_coords_rounded = precision * np.round(z_coords / precision)

    z_unique = np.unique(z_coords_rounded)

    distance = np.diff(z_unique)

    print(f"Layer Distance of {structure_fn}:")
    print(f"  {distance}")

    return distance, len(z_unique)


def layers_distance_diff(
    structure1_fn: str,
    structure2_fn: str,
    precision: float = 0.001,
):
    """
    统计原子层间距变化
    仅限表面构型弛豫前后的原子层间距变化
    """

    distance1, layer1_count = get_distance(structure1_fn, precision)
    distance2, layer2_count = get_distance(structure2_fn, precision)

    if layer1_count != layer2_count:
        raise ValueError("Layer count of two structures is different! Exit.")
    else:
        distance_diff = distance1 - distance2

        print(
            f"\nLayer distance difference between {structure1_fn} with {structure2_fn}:"
        )
        print(f"  {distance_diff}")


if __name__ == "__main__":
    structure1_fn = sys.argv[1]
    structure2_fn = sys.argv[2]

    layers_distance_diff(
        structure1_fn,
        structure2_fn,
    )
