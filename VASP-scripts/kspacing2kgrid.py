#!/usr/bin/env python3

"""
Get K-grid values from K-spacing value and POSCAR

reference: https://github.com/Tonner-Zech-Group/VASP-tools/blob/main/src/tools4vasp/kspacing2kgrid.py
"""

import argparse

import numpy as np
from ase.io import read


def get_kgrid(kspacing: float, ouput: bool = False):
    """Get K-grid values from K-spacing value and POSCAR"""

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
    print("K-grid for K-spacing {} Å^-1: {} {} {}".format(kspacing, *kgrid))

    if ouput:
        with open("KPOINTS", "w") as f:
            f.write(f"K-Spacing Value to Generate K-Mesh: {kspacing:.2f}\n")
            f.write("0\n")
            f.write("Gamma\n")
            f.write("  ".join(map(str, kgrid)))
            f.write("\n")
            f.write("0.0  0.0  0.0\n")

        print("\nKPOINTS generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get KGRID values from KSPACING value and POSCAR file.",
    )

    parser.add_argument("kspacing", type=float, help="K-spacing value")
    parser.add_argument("-o", "--output", action="store_true", help="Generate KPOINTS")

    args = parser.parse_args()

    get_kgrid(args.kspacing, args.output)
