#!/usr/bin/env python3

"""
获取结构的 Pearson 符号

reference: https://www.ctcms.nist.gov/potentials/iprPy/notebook/crystal_space_group.html
"""

import argparse

import spglib
from ase.atoms import Atoms
from ase.io import read


def get_pearson_symbol(structure_fn: str = "POSCAR"):
    """获取结构的 Pearson 符号"""

    atoms: Atoms = read(structure_fn)
    natoms = len(atoms)

    lattice = atoms.cell
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell = (lattice, positions, numbers)

    sym_data = spglib.get_symmetry_dataset(cell)

    spg_type = spglib.get_spacegroup_type(sym_data.hall_number)

    # 晶系
    if spg_type.number <= 2:
        crystal_system = "a"
    elif spg_type.number <= 15:
        crystal_system = "m"
    elif spg_type.number <= 74:
        crystal_system = "o"
    elif spg_type.number <= 142:
        crystal_system = "t"
    elif spg_type.number <= 194:
        crystal_system = "h"
    else:
        crystal_system = "c"

    # 布拉维点阵类型；若为 A/B 面的有心化，转为 C 面
    latticetype = spg_type.international[0]
    if latticetype in ["A", "B"]:
        latticetype = "C"

    pearson_symbol = crystal_system + latticetype + str(natoms)

    return pearson_symbol


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get Pearson symbol of a structure",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn", nargs="?", default="POSCAR", help="structure filename"
    )

    args = parser.parse_args()

    pearson_symbol = get_pearson_symbol(structure_fn=args.structure_fn)
    print(f"Pearson Symbol: {pearson_symbol}")
