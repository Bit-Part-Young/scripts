#!/usr/bin/env python3

"""识别原子层及其对应的原子"""

import argparse

import numpy as np
import pandas as pd
from ase.io import read


def layer_identify(
    structure_fn: str,
    layer_indices: int | list[int] | None = None,
    precision: float = 0.001,
    axis: str = "z",
):
    """识别原子层及其对应的原子"""

    if structure_fn.split(".")[-1] in ["lmp", "data", "lammps-data"]:
        atoms = read(structure_fn, format="lammps-data", sort_by_id=False)
    elif structure_fn.split(".")[-1] in ["xyz", "extxyz"]:
        atoms = read(structure_fn, format="extxyz")
    else:
        atoms = read(structure_fn)

    element_list = atoms.get_chemical_symbols()

    positions = atoms.get_positions(wrap=True)
    positions_frac = atoms.get_scaled_positions(wrap=True)

    if axis == "z":
        positions_z = positions[:, 2]
    elif axis == "y":
        positions_z = positions[:, 1]
    elif axis == "x":
        positions_z = positions[:, 0]
    else:
        raise ValueError(f"Invalid axis: {axis}")

    positions_z_rounded = precision * np.round(positions_z / precision)
    positions_z_unique = np.unique(positions_z_rounded)

    print(f"Atoms count: {len(atoms)}; Layers count: {len(positions_z_unique)}.\n")

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
        "structure_fn", nargs="?", default="POSCAR", help="structure filename"
    )

    parser.add_argument(
        "-li",
        "--layer_indices",
        type=int,
        nargs="*",
        metavar="N",
        help="atomic layer indices",
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
    layer_indices = args.layer_indices
    axis = args.axis

    layer_identify(structure_fn=structure_fn, layer_indices=layer_indices, axis=axis)
