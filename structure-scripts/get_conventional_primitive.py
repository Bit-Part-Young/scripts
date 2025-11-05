#!/usr/bin/env python3

"""获取给定构型的原胞 primitive、单胞 conventional"""

import argparse

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def get_conventional_primitive(structure_fn: str, mode: str, prec: float):

    structure_initial = Poscar.from_file(structure_fn).structure
    ss = SpacegroupAnalyzer(structure_initial, symprec=prec)

    if mode in ["primitive", "prim"]:
        structure_final = ss.get_primitive_standard_structure()
        output_fn = "primitive.vasp"

        print("\nThe primitive structure:\n")
    elif mode in ["conventional", "conv"]:
        structure_final = ss.get_conventional_standard_structure()
        output_fn = "conventional.vasp"

        print("\nThe conventional structure:\n")

    poscar = Poscar(structure_final)
    poscar.write_file(output_fn, significant_figures=12)

    print(structure_final)

    print(f"\n{output_fn} is generated.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get the conventional/primitive structure.",
        allow_abbrev=True,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "mode",
        choices=["primitive", "prim", "conventional", "conv"],
        help="To define the mode of the structure; by default: primitive.",
    )
    parser.add_argument(
        "-i",
        "--input",
        default="POSCAR",
        metavar="FILE",
        help="input structure filename",
    )

    parser.add_argument(
        "-p",
        "--prec",
        type=float,
        default=0.001,
        metavar="FLOAT",
        help="symmetry precision",
    )

    args = parser.parse_args()

    get_conventional_primitive(mode=args.mode, structure_fn=args.input, prec=args.prec)
