#!/usr/bin/env python3

"""
识别每个原子所在的原子层

Usage:
layer_index_identification.py POSCAR      # 显示所有原子的原子层
layer_index_identification.py POSCAR 3    # 显示第 3 个原子层为的所有原子
layer_index_identification.py POSCAR 3 4  # 显示第 3、4 个原子层为的所有原子
"""

import sys

import numpy as np
import pandas as pd
from pymatgen.core.structure import Structure


def layer_index_identification(
    structure_fn: str,
    layer_number: int | list[int] | None = None,
    precision: float = 0.001,
):
    """识别每个原子所在的原子层"""

    structure = Structure.from_file(structure_fn)

    positions = structure.cart_coords

    z_coords = positions[:, 2]
    z_coords_rounded = precision * np.round(z_coords / precision)

    z_unique = np.unique(z_coords_rounded)
    z_unique_sorted = np.sort(z_unique)

    print(f"Number of atom layer: {len(z_unique)}.\n")

    # 获取每个原子所在的原子层
    layer_num_list = []
    for z_coord in z_coords_rounded:
        for i, z in enumerate(z_unique_sorted, start=1):
            if np.isclose(z, z_coord):
                layer_num_list.append(i)

    data = {
        "x": positions[:, 0],
        "y": positions[:, 1],
        "z": positions[:, 2],
        "layer_number": layer_num_list,
    }
    df = pd.DataFrame(data)

    if layer_number is None:
        print(df)
    elif isinstance(layer_number, int):
        print(df[df["layer_number"] == layer_number])
    elif isinstance(layer_number, list):
        print(df[df["layer_number"].isin(layer_number)])


if __name__ == "__main__":

    structure_fn = sys.argv[1]
    layer_number = int(sys.argv[2]) if len(sys.argv) == 3 else None

    if len(sys.argv) > 3:
        layer_number = [int(arg) for arg in sys.argv[2:]]

    layer_index_identification(
        structure_fn=structure_fn,
        layer_number=layer_number,
    )
