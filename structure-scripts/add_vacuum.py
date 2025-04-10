#!/usr/bin/env python3

"""添加真空层（包括：z 方向顶部、底部、顶部底部两端）"""

import argparse

from ase.atoms import Atoms
from ase.io import read, write


def add_vacuum(atoms: Atoms, vacuum: float, mode: str) -> Atoms:
    """在 z 方向顶部、底部、顶部底部两端添加真空层"""

    cell = atoms.get_cell()
    positions = atoms.get_positions()

    if mode == "top":
        cell[2, 2] += vacuum
    elif mode == "bottom":
        cell[2, 2] += vacuum
        positions[:, 2] += vacuum
        atoms.set_positions(positions)
    elif mode == "both":
        cell[2, 2] += vacuum * 2
        positions[:, 2] += vacuum
        atoms.set_positions(positions)

    atoms.set_cell(cell)

    return atoms


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Add vacuum to the top, bottom, or both ends in z axis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument("mode", choices=["top", "bottom", "both"], help="Mode to add vacuum")
    parser.add_argument("vacuum", type=float, default=15.0, help="Vacuum size")
    parser.add_argument("structure_fn", type=str, default="POSCAR", help="Structure filename")

    args = parser.parse_args()

    atoms = read(args.structure_fn)
    atoms = add_vacuum(atoms, args.vacuum, args.mode)

    write("add_vacuum.vasp", atoms, format="vasp", vasp5=True, direct=True)
    print(f"Added {args.vacuum} Å vacuum to the {args.mode} and saved to add_vacuum.vasp.")
