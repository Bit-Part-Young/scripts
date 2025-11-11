#!/usr/bin/env python3

"""
获取结构的晶体学对称性信息

已知问题: 对于 Ti3P 型 Nb3Si 结构，其 wyckoff 为 8g NbI、8g NbII、8g NbIII，程序将其统计成了 24g，是错误的
解决方式: 添加 -v/--verbose 选项，进一步查看 Equivalent Site Group 信息进行检查

reference: https://github.com/nanyanshouhu/Defect_generator/blob/main/wyckoff_position_finder.py
"""

import argparse
import os
from collections import Counter

import spglib
from ase.atoms import Atoms
from ase.io import read
from pymatgen.core import Structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def get_pearson_symbol(atoms: Atoms):
    """获取结构的 Pearson 符号"""

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


def symmetry_info(
    structure_fn: str = "POSCAR",
    symprec: float = 0.01,
    angle_tolerance: float = 5,
    verbose: bool = False,
):
    """获取结构的晶体学对称性信息"""

    input_fn_basename = os.path.basename(structure_fn)
    input_format = input_fn_basename.split(".")[-1]

    if input_format in ["POSCAR", "CONTCAR", "vasp"]:
        structure = Structure.from_file(structure_fn)
        atoms = structure.to_ase_atoms()
    elif input_format in ["xyz", "extxyz"]:
        atoms = read(structure_fn, format="extxyz")
        structure = Structure.from_ase_atoms(atoms)
    elif input_format in ["lmp", "data", "lammps-data"]:
        # sort_by_id 为 True 时，原子的 id 不连续时，会出现报错；为 False，不会报错，建议设为 False
        atoms = read(structure_fn, format="lammps-data", sort_by_id=False)
        structure = Structure.from_ase_atoms(atoms)

    pearson_symbol = get_pearson_symbol(atoms)

    sga = SpacegroupAnalyzer(
        structure=structure, symprec=symprec, angle_tolerance=angle_tolerance
    )

    symmetry_dataset = sga.get_symmetry_dataset()
    wyckoff_positions = symmetry_dataset.wyckoffs
    space_group_number = symmetry_dataset.number
    space_group_symbol = symmetry_dataset.international

    print(f"Space Group Symbol: {space_group_symbol}")
    print(f"Space Group Number: {space_group_number}")
    print(f"Pearson Symbol: {pearson_symbol}")

    # 统计每种元素的 Wyckoff positions
    wyckoff_count = {}
    for site, wyckoff in zip(structure, wyckoff_positions):
        site: PeriodicSite
        element = site.species_string
        if element not in wyckoff_count:
            wyckoff_count[element] = []
        wyckoff_count[element].append(wyckoff)

    # Wyckoff 信息
    print("\nWyckoff info:")
    for element, positions in wyckoff_count.items():
        position_count = Counter(positions)
        wyckoff_summary = " ".join(
            [f"{count}{wyckoff}" for wyckoff, count in position_count.items()]
        )
        print(f"{element}: {wyckoff_summary}")

    if verbose:
        # 获取对称性等同的原子位置
        equivalent_sites = sga.get_symmetrized_structure().equivalent_sites
        print("\nEquivalent Sites:")
        for i, sites in enumerate(equivalent_sites):
            print(f"\nEquivalent Site Group {i + 1}:")
            for site in sites:
                print(f"  - {site.species_string} at {site.frac_coords.round(5)}")
    else:
        print("\nWyckoff info maybe incorrect.", end=" ")
        print("Use -v/--verbose option to check Equivalent Site Group info.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get crystal symmetry info of a structure.", epilog="Author: SLY."
    )

    parser.add_argument(
        "structure_fn", nargs="?", default="POSCAR", help="structure filename"
    )

    poscar_fn = parser.add_argument(
        "--symprec", type=float, default=0.01, help="symmetry precision"
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="show equivalent site group info"
    )

    args = parser.parse_args()

    symmetry_info(
        structure_fn=args.structure_fn, symprec=args.symprec, verbose=args.verbose
    )
