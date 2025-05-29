#!/usr/bin/env python3

"""给 model.xyz 添加 group 信息；group method 原子 z 轴坐标值"""

import argparse

import numpy as np
from ase.io import read, write


def add_groups(
    structure_fn: str,
    output_fn: str = "model.xyz",
):

    if structure_fn.split(".")[-1] == "xyz":
        atoms = read(structure_fn, format="extxyz")
    else:
        atoms = read(structure_fn)

    group_list: list[int] = []
    for atom in atoms:
        if atom.position[2] < atoms.cell[2, 2] / 2:
            group_list.append(0)
        else:
            group_list.append(1)

    group_array = np.array(group_list)
    atoms.new_array("group", group_array)

    write(output_fn, atoms, format="extxyz")

    print(f"\nGroup info added to {structure_fn} and saved to {output_fn}!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add groups info to structure according to atom z-axis coordinate.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        help="structure filename",
    )

    parser.add_argument(
        "-o",
        "--output_fn",
        type=str,
        default="model.xyz",
        help="output filename",
    )

    args = parser.parse_args()
    structure_fn = args.structure_fn
    output_fn = args.output_fn

    add_groups(structure_fn, output_fn)
