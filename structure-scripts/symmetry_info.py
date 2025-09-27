#!/usr/bin/env python3

"""
获取结构的对称性信息

已知问题: 对于 Ti3P 型 Nb3Si 结构，其 wyckoff 为 8g NbI、8g NbII、8g NbIII，程序将其统计成了 24g，是错误的
解决方式: 添加 -v/--verbose 选项，进一步查看 Equivalent Site Group 信息进行检查

reference: https://github.com/nanyanshouhu/Defect_generator/blob/main/wyckoff_position_finder.py
"""

import argparse
from collections import Counter

from ase.io import read
from ase.io.vasp import read_vasp
from pymatgen.core import Structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def symmetry_info(
    structure_fn: str = "POSCAR",
    symprec: float = 0.01,
    angle_tolerance: float = 5,
    verbose: bool = False,
):
    """获取结构的对称性信息"""

    structure = Structure.from_file(structure_fn)

    sga = SpacegroupAnalyzer(
        structure=structure, symprec=symprec, angle_tolerance=angle_tolerance
    )

    symmetry_dataset = sga.get_symmetry_dataset()
    wyckoff_positions = symmetry_dataset.wyckoffs
    space_group_number = symmetry_dataset.number
    space_group_symbol = symmetry_dataset.international

    print(f"Space Group Symbol: {space_group_symbol}")
    print(f"Space Group Number: {space_group_number}")

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
        description="Get structure symmetry info.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn", nargs="?", default="POSCAR", help="structure filename"
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="show equivalent site group info"
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn
    symmetry_info(structure_fn=structure_fn, verbose=args.verbose)
