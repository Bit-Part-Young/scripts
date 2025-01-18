#!/usr/bin/env python3

"""
Get KSPACING values from KPOINTS and POSCAR files

reference: https://github.com/Tonner-Zech-Group/VASP-tools/blob/main/src/tools4vasp/kspacing2kgrid.py
"""

import argparse

import numpy as np
from ase.io import read


def get_kspacing():
    """Get KSPACING values from KPOINTS and POSCAR files"""

    atoms = read("POSCAR")
    cellparams = np.array(atoms.cell.cellpar()).round(3)

    print("Cell parameters: a={} b={} c={} alpha={} beta={} gamma={}".format(*cellparams))

    # 不包含 2*pi 因子
    cell_reciprocal = atoms.cell.reciprocal()[:].round(5)
    print("Reciprocal cell (without 2*pi):")
    print(cell_reciprocal)

    with open("KPOINTS", "r") as f:
        kpoints = f.readlines()

    if kpoints[1].strip() != "0":
        raise ValueError("KPOINTS file does not use Automatic Scheme, but only this is supported!")

    kgrid = [int(k) for k in kpoints[3].split()]
    assert len(kgrid) == 3, "Expected to find 3 integers in line 4 of KPOINTS file!"
    print("KGRID in KPOINTS: {} {} {}".format(*kgrid))

    kspacing = [np.linalg.norm(cell_reciprocal[i]) * 2 * np.pi / kgrid[i] for i in range(3)]
    kspacing = np.array(kspacing).round(3)
    print("KSPACING: {} Å^-1, {} Å^-1, {} Å^-1".format(*kspacing))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get KSPACING values from POSCAR and KPOINTS files.",
    )

    args = parser.parse_args()

    get_kspacing()
