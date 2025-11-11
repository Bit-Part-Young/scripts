#!/usr/bin/env python3

"""
生成 Ga2Hf (r-Al2Zr) 结构
Date: 2025-11-10

晶体学信息 reference:
- First-principles calculation of structural energetics of Al-TM (TM=Ti, Zr, Hf) intermetallics
"""

from ase.atoms import Atoms
from ase.io import write
from ase.spacegroup import crystal
from pymatgen.core.structure import Structure


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


def structure_r_Al2Ti_generation():
    """生成 C14 结构"""

    a, b, c = 3.97, 3.97, 24.32
    spacegroup = 141

    wyckoff_positions = [
        ("Al", [0.0, 0.0, 0.25]),
        ("Al", [0.00, 0.00, 0.4100]),
        ("Ti", [0.00, 0.00, 0.0770]),
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

    print("\nr-Al2Ti structure generated.")


if __name__ == "__main__":
    structure_r_Al2Ti_generation()
