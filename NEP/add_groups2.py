#!/usr/bin/env python3

"""给 model.xyz 添加 group 信息；group method 为 元素种类顺序"""

import argparse

import numpy as np
from ase.io import read, write


def add_groups(
    structure_fn: str,
    elements: list[str],
    output_fn: str = "model.xyz",
):

    if structure_fn.split(".")[-1] == "xyz":
        atoms = read(structure_fn, format="extxyz")
    else:
        atoms = read(structure_fn)

    group_list: list[int] = []
    for element in atoms.get_chemical_symbols():
        if element in elements:
            group = elements.index(element)
        else:
            raise ValueError(
                f"Element {element} not found in the provided elements list"
            )
        group_list.append(group)

    group_array = np.array(group_list)
    atoms.new_array("group", group_array)

    write(output_fn, atoms, format="extxyz")

    print(f"\nGroup info added to {structure_fn} and saved to {output_fn}!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add groups info to structure according to elements sequence.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        help="structure filename",
    )
    parser.add_argument(
        "-ess",
        type=str,
        nargs="+",
        required=True,
        help="element symbol sequence, e.g. Ti Al Nb",
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
    elements = args.ess
    output_fn = args.output_fn

    add_groups(structure_fn, elements, output_fn)
