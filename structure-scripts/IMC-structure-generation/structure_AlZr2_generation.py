#!/usr/bin/env python3

"""
生成 hP6 B8_2 (AlZr2) 结构
Date: 2025-12-03

晶体学信息 reference:
- First-principles calculation of structural energetics of Al-TM (TM=Ti, Zr, Hf) intermetallics
"""

import numpy as np
from ase.atoms import Atoms
from ase.io import write
from ase.spacegroup import crystal
from pymatgen.core.structure import Structure

np.set_printoptions(precision=5, suppress=True)


def get_structure_info(structure: Structure):
    """输出原子位点信息"""

    print(f"\nnatoms: {len(structure)}.")
    print("\natoms site info:")

    for site in structure:
        print(site.species, site.frac_coords)


def get_atoms_info(atoms: Atoms):
    """输出原子位点信息"""

    print(f"\nnatoms: {len(atoms)}.")
    print("\natoms site info:")

    for atom in atoms:
        print(atom.symbol, atom.scaled_position)


def structure_AlZr2_generation():
    """生成 hP6 (AlZr2) 结构"""

    a, b, c = 4.890, 4.890, 5.900
    spacegroup = 194

    wyckoff_positions = [
        ("Al", [0.3333, 0.6667, 0.25]),
        ("Zr", [0.0, 0.0, 0.5]),
        ("Zr", [0.3333, 0.6667, 0.75]),
    ]

    atoms = crystal(
        symbols=[wyckoff_position[0] for wyckoff_position in wyckoff_positions],
        basis=[wyckoff_position[1] for wyckoff_position in wyckoff_positions],
        spacegroup=spacegroup,
        cellpar=[a, b, c, 90, 90, 120],
        symprec=0.01,
    )

    get_atoms_info(atoms)

    write("POSCAR", atoms, format="vasp", direct=True, sort=True)

    print("\nAl4Zr5 structure generated.")


if __name__ == "__main__":
    structure_AlZr2_generation()
