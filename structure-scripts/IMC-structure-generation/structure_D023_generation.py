#!/usr/bin/env python3

"""
生成 D023 (Al3Zr) 结构
Date: 2025-11-04

晶体学信息 reference:
- First-principles calculation of structural energetics of Al-TM (TM=Ti, Zr, Hf) intermetallics
"""

from ase.atoms import Atoms
from ase.io import write
from ase.spacegroup import crystal
from pymatgen.core.structure import Structure


def get_structure_info(structure: Structure):
    """输出原子位点信息"""

    print(f"natoms: {len(structure)}")
    print("\natoms site info:")

    for site in structure:
        print(site.species, site.frac_coords)


def get_atoms_info(atoms: Atoms):
    """输出原子位点信息"""

    print(f"\nnatoms: {len(atoms)}.")
    print("\natoms site info:")

    for atom in atoms:
        print(atom.symbol, atom.scaled_position)


def structure_D023_generation():
    """生成 D023 结构"""

    a, b, c = 4.000, 4.000, 17.300
    spacegroup = 139

    wyckoff_positions = [
        ("Al", [0.0, 0.5, 0.0]),
        ("Al", [0.0, 0.5, 0.25]),
        ("Al", [0.0, 0.0, 0.375]),
        ("Zr", [0.0, 0.0, 0.119]),
    ]

    atoms = crystal(
        symbols=[wyckoff_position[0] for wyckoff_position in wyckoff_positions],
        basis=[wyckoff_position[1] for wyckoff_position in wyckoff_positions],
        spacegroup=spacegroup,
        cellpar=[a, b, c, 90, 90, 90],
        symprec=0.01,
    )

    get_atoms_info(atoms)

    write("POSCAR", atoms, format="vasp", direct=True, sort=True)

    print("\nD023 structure generated.")


if __name__ == "__main__":
    structure_D023_generation()
