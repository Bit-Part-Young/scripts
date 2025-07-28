#!/usr/bin/env python3

"""识别原子层及其对应的原子"""

import argparse

import numpy as np
import pandas as pd
from pymatgen.core.structure import Structure


def layer_identify(
    structure_fn: str,
    layer_indices: int | list[int] | None = None,
    precision: float = 0.001,
):
    """识别原子层及其对应的原子"""

    structure = Structure.from_file(structure_fn)

    element_list = [site.species_string for site in structure]

    positions = structure.cart_coords
    positions_frac = structure.frac_coords

    positions_z = positions[:, 2]
    positions_z_rounded = precision * np.round(positions_z / precision)
    positions_z_unique = np.unique(positions_z_rounded)

    print(
        f"Atoms count: {structure.num_sites}; Layers count: {len(positions_z_unique)}.\n"
    )

    # 获取每个原子所在的原子层
    layer_index_list = []
    for z_coord in positions_z_rounded:
        for i, z_unique in enumerate(positions_z_unique, start=1):
            if np.isclose(z_unique, z_coord):
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

    if layer_indices is None:
        print(df)
    elif isinstance(layer_indices, int):
        print(df[df["layer_index"] == layer_indices])
    elif isinstance(layer_indices, list):
        print(df[df["layer_index"].isin(layer_indices)])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Identify atomic layers.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
    )

    parser.add_argument(
        "structure_fn",
        nargs="?",
        default="POSCAR",
        metavar="structure_fn",
        help="structure filename",
    )

    parser.add_argument(
        "-li",
        "--layer_indices",
        type=int,
        nargs="*",
        metavar="layer_indices",
        help="atomic layer indices",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn
    layer_indices = args.layer_indices

    layer_identify(
        structure_fn=structure_fn,
        layer_indices=layer_indices,
    )
