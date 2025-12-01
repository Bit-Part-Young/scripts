#!/usr/bin/env python3

"""
生成 NEB 插值构型（使用 ASE NEB 模块）

reference: https://github.com/DYang90/ase_extend/blob/main/makeneb.py
"""

import argparse
import os

from ase.atoms import Atoms
from ase.io import read, write
from ase.mep import NEB


def neb_interpolate(
    initial_structure_fn: str,
    final_structure_fn: str,
    climb: bool = False,
    nimages: int = 7,
    neb_method: str = "aseneb",
    interpolate_method: str = "idpp",
):

    atoms_initial = read(initial_structure_fn)
    atoms_final = read(final_structure_fn)

    images: list[Atoms] = [atoms_initial]
    images += [atoms_initial.copy() for _ in range(nimages)]
    images += [atoms_final]

    # NEB 类的 climb method 参数值应该对插值构型的生成无影响？
    if neb_method in ["spline", "string"]:
        neb = NEB(
            images,
            climb=climb,
            method=neb_method,
            allow_shared_calculator=True,
            precon="Exp",
        )
    else:
        neb = NEB(images, climb=climb, method=neb_method, allow_shared_calculator=True)

    if interpolate_method == "linear":
        neb.interpolate(method="linear", apply_constraint=False)
    elif interpolate_method == "idpp":
        neb.interpolate(method="idpp", apply_constraint=False)

    return images


def write_guess(images: list[Atoms], output: bool = False):

    for index, image in enumerate(images):
        folder = f"{index:02d}"
        os.makedirs(folder, exist_ok=True)
        write(f"{folder}/POSCAR", image, format="vasp", direct=True)

    if output:
        xyz_filename = "neb.xyz"
        if os.path.exists(xyz_filename):
            os.remove(xyz_filename)
        write(xyz_filename, images, format="extxyz", append=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=" NEB Interpolation with ASE package.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "-i",
        "--images",
        nargs=2,
        metavar="FILE",
        help="initial and final structure filenames",
    )
    parser.add_argument(
        "--climb",
        action="store_true",
        help="use climbing image",
    )
    parser.add_argument(
        "-im",
        "--interpolate_method",
        type=str,
        choices=["linear", "idpp"],
        default="idpp",
        metavar="STR",
        help="interpolate method",
    )
    parser.add_argument(
        "-ni",
        "--nimages",
        type=int,
        default=7,
        metavar="N",
        help="number of interpolation images",
    )

    args = parser.parse_args()

    images = neb_interpolate(
        initial_structure_fn=args.images[0],
        final_structure_fn=args.images[1],
        climb=args.climb,
        nimages=args.nimages,
        interpolate_method=args.interpolate_method,
    )

    write_guess(images, output=True)
