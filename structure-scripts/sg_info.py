#!/usr/bin/env python3

"""
获取结构的空间群信息

reference: https://github.com/ajjackson/mctools/blob/develop/mctools/generic/get_spacegroup.py
"""

import argparse

from ase.io import read
from spglib import get_spacegroup


def sg_info(structure_fn: str = "POSCAR"):
    """获取结构的空间群信息"""

    atoms = read(structure_fn)
    cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)

    print("| Threshold / Å |    Space group    |")
    print("|---------------|-------------------|")

    threshold_list = [1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]

    for threshold in threshold_list:
        print(
            "|    {0:0.5f}    |  {1: ^16} |".format(
                threshold, get_spacegroup(cell, symprec=threshold)
            )
        )

    print("|---------------|-------------------|")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get space group info.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        nargs="?",
        default="POSCAR",
        help="Structure filename.",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn

    sg_info(structure_fn=structure_fn)
