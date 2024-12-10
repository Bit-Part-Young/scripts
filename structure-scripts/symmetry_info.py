"""
获取结构的对称性信息

reference: https://github.com/nanyanshouhu/Defect_generator/blob/main/wyckoff_position_finder.py
"""

import argparse
from collections import Counter

from ase.io import read
from pymatgen.core import Structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def symmetry_info(
    structure_fn: str = "POSCAR",
    symprec: float = 0.01,
    angle_tolerance: float = 5,
):
    """获取结构的对称性信息"""

    ase_formats = ["xyz", "lammps-data", "lmp"]

    structure_fn_format = structure_fn.split(".")[-1]
    if structure_fn_format in ase_formats:
        if structure_fn_format == "lmp":
            atoms = read(structure_fn, format="lammps-data")
        else:
            atoms = read(structure_fn, format=structure_fn_format)

        structure = Structure.from_ase_atoms(atoms)
    else:
        structure = Structure.from_file(structure_fn)

    sga = SpacegroupAnalyzer(
        structure=structure,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get structure symmetry info.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        nargs="?",
        default="POSCAR",
        help="Structure filename",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn
    symmetry_info(structure_fn=structure_fn)
