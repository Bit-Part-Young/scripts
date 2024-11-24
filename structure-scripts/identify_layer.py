#!/usr/bin/env python3

"""
识别每个原子所在的原子层

Usage:
    identify_layer.py POSCAR      # 显示所有原子的原子层
    identify_layer.py POSCAR 3    # 显示第 3 个原子层为的所有原子
    identify_layer.py POSCAR 3 4  # 显示第 3、4 个原子层为的所有原子
"""

import sys

import numpy as np
import pandas as pd
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure


def identify_layer(
    structure_fn: str,
    layer_index: int | list[int] | None = None,
    precision: float = 0.001,
):
    """
    识别每个原子所在的原子层
    """

    structure = Structure.from_file(structure_fn)

    element_list = structure.species
    element_list = [Element(element).symbol for element in element_list]

    positions = structure.cart_coords
    direct_positions = structure.frac_coords

    z_coords = positions[:, 2]
    z_coords_rounded = precision * np.round(z_coords / precision)

    z_unique = np.unique(z_coords_rounded)

    print(f"Layers count: {len(z_unique)}.\n")

    # 获取每个原子所在的原子层
    layer_index_list = []
    for z_coord in z_coords_rounded:
        for i, z in enumerate(z_unique, start=1):
            if np.isclose(z, z_coord):
                layer_index_list.append(i)

    data = {
        "x": positions[:, 0],
        "y": positions[:, 1],
        "z": positions[:, 2],
        "layer_index": layer_index_list,
        "element": element_list,
        "xs": direct_positions[:, 0],
        "ys": direct_positions[:, 1],
        "zs": direct_positions[:, 2],
    }
    df = pd.DataFrame(data).round(5)

    if layer_index is None:
        print(df)
    elif isinstance(layer_index, int):
        print(df[df["layer_index"] == layer_index])
    elif isinstance(layer_index, list):
        print(df[df["layer_index"].isin(layer_index)])


if __name__ == "__main__":

    structure_fn = sys.argv[1]
    layer_index = int(sys.argv[2]) if len(sys.argv) == 3 else None

    if len(sys.argv) > 3:
        layer_index = [int(arg) for arg in sys.argv[2:]]

    identify_layer(
        structure_fn=structure_fn,
        layer_index=layer_index,
    )
