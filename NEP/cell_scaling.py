#!/usr/bin/env python3

"""对结构进行 cell 缩放"""

import argparse

from ase.io import read, write


def cell_scaling(
    structure_fn: str = "POSCAR",
    scaling: float = 1.05,
    output_fn: str = "cell_scaled.vasp",
):
    """对结构进行 cell 缩放"""

    atoms = read(structure_fn, format="vasp")

    atoms_copy = atoms.copy()

    scaled_positions = atoms.get_scaled_positions()

    atoms.cell *= scaling

    atoms.set_scaled_positions(scaled_positions)

    write(output_fn, atoms, format="vasp", direct=True, sort=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Structure cell scaling.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn",
        default="POSCAR",
        help="input structure filename",
    )
    parser.add_argument(
        "scaling",
        type=float,
        default=1.05,
        metavar="FLOAT",
        help="cell scaling factor",
    )
    parser.add_argument(
        "-o",
        "--output_fn",
        default="cell_scaled.vasp",
        metavar="FILE",
        help="output structure filename",
    )

    args = parser.parse_args()

    cell_scaling(args.structure_fn, args.scaling, args.output_fn)
