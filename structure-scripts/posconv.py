#!/usr/bin/env python3

"""构型文件格式转换"""

import argparse

from ase.atoms import Atoms
from ase.io import read, write


def posconv(
    input_fn: str,
    output_fn: str,
):
    """构型文件格式转换"""

    atoms: Atoms = read(input_fn)

    print(atoms.get_chemical_symbols())

    output_format = output_fn.split(".")[-1]
    if output_format in ["vasp", "POSCAR"]:
        write(
            output_fn,
            images=atoms,
            format="vasp",
            direct=True,
            sort=True,
        )
    elif output_format in ["lammps-data", "lmp"]:
        write(
            output_fn,
            images=atoms,
            format="lammps-data",
            atom_style="atomic",
            masses=True,
        )
    else:
        write(
            output_fn,
            images=atoms,
            format=output_format,
        )

    print(f"\nConvert {input_fn} to {output_fn}!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Structure file format convert. Support most ASE recognized formats.",
        epilog="Author: SLY.",
        allow_abbrev=True,
    )

    parser.add_argument(
        "input_fn",
        type=str,
        help="Input filename, include: POSCAR, CONTCAR, *.vasp, *.xsd, *.xyz, *.lammps-data, ...",
    )

    parser.add_argument(
        "output_fn",
        type=str,
        help="Output filename, include: POSCAR, *.vasp, *.xsd, *.xyz, *.lammps-data, ...",
    )

    args = parser.parse_args()

    input_fn = args.input_fn
    output_fn = args.output_fn

    posconv(
        input_fn=input_fn,
        output_fn=output_fn,
    )
