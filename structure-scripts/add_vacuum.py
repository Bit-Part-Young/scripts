#!/usr/bin/env python3

"""在指定轴（默认 z 轴）顶部、底部、顶部底部两端添加真空层"""

import argparse

from ase.atoms import Atoms
from ase.io import read, write


def add_vacuum(
    structure_fn: str,
    vacuum: float,
    mode: str,
    axis: str = "z",
    save_poscar: bool = False,
) -> Atoms:
    """在指定轴（默认 z 轴）顶部、底部、顶部底部两端添加真空层"""

    atoms = read(structure_fn, format="vasp")

    cell = atoms.get_cell()
    positions = atoms.get_positions()

    if axis == "z":
        axis_index = 2
    elif axis == "y":
        axis_index = 1
    elif axis == "x":
        axis_index = 0

    if mode == "top":
        cell[axis_index, axis_index] += vacuum
    elif mode == "bottom":
        cell[axis_index, axis_index] += vacuum
        positions[:, axis_index] += vacuum
        atoms.set_positions(positions)
    elif mode == "both":
        cell[axis_index, axis_index] += vacuum * 2
        positions[:, axis_index] += vacuum
        atoms.set_positions(positions)

    atoms.set_cell(cell)

    if save_poscar:
        output_fn = "vacuum.vasp"
        write(output_fn, atoms, format="vasp", vasp5=True, direct=True)
        print(
            f"Added {args.vacuum} Å vacuum to {args.mode} side and saved to {output_fn}."
        )
    else:
        print(f"Added {args.vacuum} Å vacuum to {args.mode} side.")

    return atoms


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Add vacuum to the top, bottom, or both ends in z axis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument(
        "structure_fn",
        nargs="?",
        default="POSCAR",
        metavar="structure_fn",
        help="structure filename",
    )
    parser.add_argument(
        "vacuum",
        type=float,
        default=10.0,
        metavar="vacuum",
        help="vacuum thickness",
    )
    parser.add_argument(
        "mode",
        choices=["top", "bottom", "both"],
        metavar="mode",
        help="mode to add vacuum",
    )

    parser.add_argument(
        "--axis",
        choices=["x", "y", "z"],
        default="z",
        metavar="axis",
        help="axis to add vacuum",
    )

    parser.add_argument("-o", action="store_true", help="write file")

    args = parser.parse_args()

    add_vacuum(args.structure_fn, args.vacuum, args.mode, args.axis, args.o)
