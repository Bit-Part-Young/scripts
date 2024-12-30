#!/usr/bin/env python3

"""将 VASP POSCAR 分数坐标范围 wrap 在 0-1 之间"""

import argparse

from ase.io import read, write


def wrap_pos(structure_fn: str):
    """将 VASP POSCAR 分数坐标 wrap 在 0-1 之间"""

    atoms = read(structure_fn)
    atoms.wrap()

    write(
        structure_fn,
        images=atoms,
        format="vasp",
        direct=True,
        sort=False,
    )

    print("The position has been wrap to 0-1.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Wrap the atomic direct coordinates with POSCAR format to 0-1.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        nargs="?",
        type=str,
        default="POSCAR",
        help="Structure filename",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn

    wrap_pos(structure_fn=structure_fn)
