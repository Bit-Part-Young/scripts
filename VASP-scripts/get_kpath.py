#!/usr/bin/env python3

"""根据构型生成 K-Path"""

import argparse

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath


# divisions = 20 vaspkit 默认值
def get_kpath(divisions: int = 20, save: bool = False):
    """根据构型生成 K-Path"""

    structure = Structure.from_file("POSCAR").to_primitive()

    kpath = HighSymmKpath(structure)

    kpoints = Kpoints.automatic_linemode(divisions=divisions, ibz=kpath)

    print(kpoints)

    if save:
        kpoints.write_file("KPOINTS")

        print("\nKPOINTS generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get K-Path from POSCAR file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "divisions",
        nargs="?",
        const=20,
        default=20,
        type=int,
        help="number of k-points along each high symmetry line",
    )

    parser.add_argument("-s", "--save", action="store_true", help="Save KPOINTS")

    args = parser.parse_args()

    get_kpath(args.divisions, args.save)
