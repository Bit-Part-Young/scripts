#!/usr/bin/env python3

"""扩胞"""

import argparse
import os

from ase.atoms import Atoms
from ase.io import read, write


def supercell_generation(
    input_fn: str,
    dimension: tuple[int, int, int],
):
    """构型文件格式转换"""

    input_fn_basename = os.path.basename(input_fn)

    atoms: Atoms = read(input_fn_basename, format="vasp")

    supercell = atoms * dimension

    output_fn = f"POSCAR.{dimension[0]}{dimension[1]}{dimension[2]}"
    write(output_fn, images=supercell, format="vasp", direct=True, sort=True)

    print(f"\n{input_fn} make supercell to {output_fn}!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Make supercell of VASP POSCAR format.",
        allow_abbrev=True,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "input_fn",
        help="input structure filename of VASP POSCAR format",
    )

    parser.add_argument(
        "-d",
        "--dim",
        type=int,
        nargs=3,
        metavar="N",
        default=[1, 1, 1],
        help="x y z dimension",
    )

    args = parser.parse_args()

    supercell_generation(input_fn=args.input_fn, dimension=args.dim)
