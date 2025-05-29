#!/usr/bin/env python3

"""给构型添加原子类型（对 MTP 势函数有用）"""

import argparse

from ase.io import read, write


def add_atom_types(structure_fn: str, elements: list[str], output_fn: str = "data.lmp"):
    """给构型添加原子类型（对 MTP 势函数有用）"""

    if structure_fn.split(".")[-1] in ["lmp", "lammps-data"]:
        atoms = read(structure_fn, format="lammps-data")
    else:
        atoms = read(structure_fn)

    write(
        output_fn,
        images=atoms,
        format="lammps-data",
        specorder=elements,
        units="metal",
        atom_style="atomic",
        masses=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add atom types to structure according to specific elements sequence.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        help="structure filename with VASP POSCAR or LAMMPS data format",
    )
    parser.add_argument(
        "-ess",
        type=str,
        nargs="+",
        required=True,
        help="element symbol sequence, e.g. Ti Al Nb",
    )
    parser.add_argument(
        "-o",
        "--output_fn",
        type=str,
        default="data.lmp",
        const="data.lmp",
        nargs="?",
        help="output filename with LAMMPS data format",
    )
    args = parser.parse_args()

    structure_fn = args.structure_fn
    elements = args.ess
    output_fn = args.output_fn

    add_atom_types(structure_fn, elements, output_fn)
