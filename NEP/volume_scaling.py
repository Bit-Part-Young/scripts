#!/usr/bin/env python3

"""生成体积缩放的结构"""

import argparse

from ase.io import read, write


def volume_scaling(
    structure_fn: str = "POSCAR",
    scaling: float = 1.05,
    output_fn: str = "volume_scaled.vasp",
):
    """生成体积缩放的结构"""

    atoms = read(structure_fn, format="vasp")

    atoms_copy = atoms.copy()

    scaled_positions = atoms.get_scaled_positions()

    atoms.cell *= scaling ** (1 / 3)

    atoms.set_scaled_positions(scaled_positions)

    write(output_fn, atoms, format="vasp", direct=True, sort=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Structure generation with cell scaling.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        default="POSCAR",
        help="structure filename (default: POSCAR)",
    )
    parser.add_argument(
        "scaling",
        type=float,
        default=1.05,
        help="cell scaling factor (default: 1.05)",
    )
    parser.add_argument(
        "-o",
        "--output_fn",
        type=str,
        default="volume_scaled.vasp",
        help="output filename",
    )

    args = parser.parse_args()

    volume_scaling(args.structure_fn, args.scaling, args.output_fn)
