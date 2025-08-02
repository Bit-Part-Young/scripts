#!/usr/bin/env python3

"""将 VASP POSCAR 中的原子分数坐标范围 wrap 在 0-1 之间"""

import argparse

from ase.io import read, write


def wrap_pos(structure_fn: str):
    """将 VASP POSCAR 中的原子分数坐标 wrap 在 0-1 之间"""

    atoms = read(structure_fn, format="vasp")
    atoms.wrap()

    write(structure_fn, images=atoms, format="vasp", direct=True, sort=False)

    print("The position has been wrap to 0-1.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Wrap the atomic direct coordinates with POSCAR format to 0-1.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn",
        nargs="?",
        type=str,
        default="POSCAR",
        metavar="structure_fn",
        help="structure filename",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn

    wrap_pos(structure_fn=structure_fn)
