#!/usr/bin/env python3

"""
获取结构的空间群信息

reference: https://github.com/ajjackson/mctools/blob/develop/mctools/generic/get_spacegroup.py
"""

import argparse
import os

from ase.io import read
from spglib import get_spacegroup


def sg_info(structure_fn: str = "POSCAR"):
    """获取结构的空间群信息"""

    input_fn_basename = os.path.basename(structure_fn)
    input_format = input_fn_basename.split(".")[-1]

    if input_format in ["POSCAR", "CONTCAR", "vasp", "PPOSCAR"]:
        atoms = read(structure_fn, format="vasp")
    elif input_format in ["xyz", "extxyz"]:
        atoms = read(structure_fn, format="extxyz")
    elif input_format in ["lmp", "data", "lammps-data"]:
        # sort_by_id 为 True 时，原子的 id 不连续时，会出现报错；为 False，不会报错，建议设为 False
        atoms = read(structure_fn, format="lammps-data", sort_by_id=False)
    elif "poscar" in structure_fn.lower():
        atoms = read(structure_fn, format="vasp")

    cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)

    print("| Threshold / Å  |    Space group    |")
    print("|----------------|-------------------|")

    threshold_list = [1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]

    for threshold in threshold_list:
        print(
            "|    {0:0.6f}    |  {1: ^16} |".format(
                threshold, get_spacegroup(cell, symprec=threshold)
            )
        )

    print("|----------------|-------------------|")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get space group info with different symprec.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn", nargs="?", default="POSCAR", help="structure filename"
    )

    args = parser.parse_args()

    sg_info(structure_fn=args.structure_fn)
