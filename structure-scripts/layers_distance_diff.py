#!/usr/bin/env python3

"""
统计原子层间距变化
仅限表面构型弛豫前后的原子层间距变化

Usage:
    layers_distance_diff.py POSCAR1 POSCAR2
"""

import sys

import numpy as np
import pandas as pd
from ase.io import read


def get_distance(
    structure_fn: str,
    precision: float = 0.001,
) -> tuple[np.ndarray, int]:
    """
    统计原子层间距
    """

    atoms = read(structure_fn)

    positions = atoms.positions

    z_coords = positions[:, 2]
    z_coords_rounded = precision * np.round(z_coords / precision)

    z_unique = np.unique(z_coords_rounded)

    distance = np.diff(z_unique)

    return distance, len(z_unique)


def layers_distance_diff(
    structure1_fn: str,
    structure2_fn: str,
    precision: float = 0.001,
) -> pd.DataFrame:
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

        data = {
            "Layer_Distance_1": distance1,
            "Layer_Distance_2": distance2,
            "Diff": distance_diff,
        }
        index = [f"{i}-{i+1}" for i in range(1, len(distance_diff) + 1)]

        df = pd.DataFrame(
            data=data,
            index=index,
        )

    return df


if __name__ == "__main__":
    structure1_fn = sys.argv[1]
    structure2_fn = sys.argv[2]

    layer_df = layers_distance_diff(
        structure1_fn,
        structure2_fn,
    )
    print(layer_df)
