#!/usr/bin/env python3

"""
生成 C15 V2Zr 结构
该方式得到的构型与文献中的相同，而 struct_gen 的则不太一致
"""

import argparse

from ase.io import write
from ase.spacegroup import crystal
from pymatgen.core.structure import Structure


# 输出原子位点信息
def get_structure_info(structure: Structure):

    print(f"natoms: {len(structure)}")
    print("\natoms site info:")

    for site in structure:
        print(site.species, site.frac_coords)


def structure_C15_generation(lattice_constant: float, element_symbols: list):
    a, b, c = lattice_constant, lattice_constant, lattice_constant
    spacegroup = 227

    wyckoff_positions = [
        (element_symbols[0], [0.0, 0.0, 0.0]),
        (element_symbols[1], [0.625, 0.625, 0.625]),
    ]

    crystal_structure = crystal(
        symbols=[wyckoff_position[0] for wyckoff_position in wyckoff_positions],
        basis=[wyckoff_position[1] for wyckoff_position in wyckoff_positions],
        spacegroup=spacegroup,
        cellpar=[a, b, c, 90, 90, 90],
        symprec=0.01,
    )

    structure = Structure.from_ase_atoms(crystal_structure)
    get_structure_info(structure)

    write("POSCAR", crystal_structure, format="vasp", direct=True, sort=True)

    print("\nC15 structure generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate C15 structure")
    parser.add_argument(
        "-e", "--element_symbols", nargs=2, metavar="STR", help="element symbols"
    )
    parser.add_argument(
        "-lc",
        "--lattice_constant",
        type=float,
        metavar="FLOAT",
        help="lattice constant",
    )
    args = parser.parse_args()

    structure_C15_generation(args.lattice_constant, args.element_symbols)
