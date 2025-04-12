#!/usr/bin/env python3

"""识别原子层及其对应的原子"""

import argparse

import numpy as np
import pandas as pd
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure


def identify_layer(
    structure_fn: str,
    layer_index: int | list[int] | None = None,
    precision: float = 0.001,
):
    """识别原子层及其对应的原子"""

    structure = Structure.from_file(structure_fn)

    element_list = structure.species
    element_list = [Element(element).symbol for element in element_list]

    positions = structure.cart_coords
    positions_frac = structure.frac_coords

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
        "xs": positions_frac[:, 0],
        "ys": positions_frac[:, 1],
        "zs": positions_frac[:, 2],
    }
    df = pd.DataFrame(data).round(5)

    if layer_index is None:
        print(df)
    elif isinstance(layer_index, int):
        print(df[df["layer_index"] == layer_index])
    elif isinstance(layer_index, list):
        print(df[df["layer_index"].isin(layer_index)])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Identify atomic layer.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
    )

    parser.add_argument(
        "structure_fn",
        nargs="?",
        type=str,
        default="POSCAR",
        help="Structure filename.",
    )

    parser.add_argument(
        "-li",
        "--layer_index",
        type=int,
        nargs="*",
        help="Atomic layer index.",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn
    layer_index = args.layer_index

    identify_layer(
        structure_fn=structure_fn,
        layer_index=layer_index,
    )
