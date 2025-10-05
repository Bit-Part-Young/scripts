#!/usr/bin/env python3

"""给定初始结构，生成其 EOS 构型，共 21 个构型，晶胞体积/晶胞参数 ±10% 范围内"""

import os

import argparse

import numpy as np
from ase.io import read, write
from typing import Literal


def eos_generation(
    structure_fn: str = "POSCAR",
    output_fn: str = "eos.xyz",
    type: Literal["volume", "cell"] = "volume",
):

    if os.path.exists(output_fn):
        os.remove(output_fn)

    atoms = read(structure_fn, format="vasp")

    scaled_positions = atoms.get_scaled_positions()

    strain_range = np.linspace(0.9, 1.1, 21)

    for strain in strain_range:
        atoms_copy = atoms.copy()

        if type == "volume":
            atoms_copy.cell *= strain ** (1 / 3)
        elif type == "cell":
            atoms_copy.cell *= strain

        atoms_copy.set_scaled_positions(scaled_positions)

        write(output_fn, atoms_copy, format="extxyz", append=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate EOS configurations of a given structure with volume variation ±10%.",
        epilog="Author: SLY.",
    )

    parser.add_argument("structure_fn", default="POSCAR", help="structure filename")
    parser.add_argument("output_fn", default="eos.xyz", help="output filename")
    parser.add_argument(
        "-t",
        "--type",
        default="volume",
        choices=["volume", "cell"],
        metavar="STR",
        help="type of strain (volume or cell)",
    )

    args = parser.parse_args()

    eos_generation(args.structure_fn, args.output_fn, args.type)
