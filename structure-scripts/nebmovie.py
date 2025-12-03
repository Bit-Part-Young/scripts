#!/usr/bin/env python3

"""生成 NEB 轨迹/movie extxyz 文件，同 VTST 中的 nebmovie.pl"""

import argparse
import os

from ase.io import read, write


def nebmovie(type: int):

    if type == 0:
        structure_fn = "POSCAR"
    elif type == 1:
        structure_fn = "CONTCAR"
    else:
        raise ValueError(f"Invalid type: {type}")

    folder_list = []
    for i in range(21):
        folder = f"{i:02d}"
        if not os.path.exists(folder):
            break
        folder_list.append(folder)

    images = []
    for index, folder in enumerate(folder_list):
        if type == 1:
            if index == 0 or index == len(folder_list) - 1:
                structure_fn = "POSCAR"
            else:
                structure_fn = "CONTCAR"

        atoms = read(f"{folder}/{structure_fn}", format="vasp")
        images.append(atoms)

    if type == 1:
        structure_fn = "CONTCAR"
    xyz_fn = f"movie_{structure_fn}.xyz"
    if os.path.exists(xyz_fn):
        os.remove(xyz_fn)
    write(xyz_fn, images, format="extxyz", append=True)

    print(f"\n{xyz_fn} generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate NEB trajectory/movie extxyz file (same as VTST's nebmovie.pl).",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "type",
        type=int,
        default=0,
        choices=[0, 1],
        help="0: POSCAR, 1: CONTCAR",
    )

    args = parser.parse_args()

    nebmovie(args.type)
