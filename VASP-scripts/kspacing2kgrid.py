#!/usr/bin/env python3

"""
Get KGRID values from KSPACING value and POSCAR file

reference: https://github.com/Tonner-Zech-Group/VASP-tools/blob/main/src/tools4vasp/kspacing2kgrid.py
"""

import argparse

import numpy as np
from ase.io import read


def get_kgrid(kspacing: float):
    """Get KGRID values from KSPACING value and POSCAR file"""

    atoms = read("POSCAR")

    cellparams = np.array(atoms.cell.cellpar()).round(3)

    print("Cell parameters: a={} b={} c={} alpha={} beta={} gamma={}".format(*cellparams))

    # 不包含 2*pi 因子
    cell_reciprocal = atoms.cell.reciprocal()[:].round(5)
    print("Reciprocal cell (without 2*pi):")
    print(cell_reciprocal)

    kgrid = [
        int(max(1, np.ceil(np.linalg.norm(cell_reciprocal[i]) * 2 * np.pi / kspacing)))
        for i in range(3)
    ]
    print("KGRID for KSPACING {} Å^-1: {} {} {}".format(kspacing, *kgrid))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get KGRID values from KSPACING value and POSCAR file.",
    )

    parser.add_argument(
        "kspacing",
        type=float,
        help="KSPACING value",
    )

    args = parser.parse_args()

    get_kgrid(args.kspacing)
