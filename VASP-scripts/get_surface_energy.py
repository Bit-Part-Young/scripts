#!/usr/bin/env python3

"""计算表面能"""

import argparse
import os

import numpy as np
from ase.io import read


def get_surface_energy(epa: float, path: str = "."):
    """计算表面能"""

    structure_fn = os.path.join(path, "POSCAR")
    outcar_fn = os.path.join(path, "OUTCAR")

    atoms = read(structure_fn)
    natoms = len(atoms)
    a, b, c = atoms.cell.lengths()
    alpha, beta, gamma = atoms.cell.angles()

    surface_area = a * b * np.sin(np.deg2rad(gamma))
    coeff = 1.6021766208 * 10

    energy = read(outcar_fn).get_potential_energy()

    surface_energy1 = (energy - epa * natoms) / (2 * surface_area)
    surface_energy2 = surface_energy1 * coeff

    print(
        f"\nsurface energy: {round(surface_energy1, 3)} eV/Å^2, {round(surface_energy2, 3)} J/m^2."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate surface energy.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "-p",
        "--path",
        default=".",
        nargs="?",
        const=".",
        metavar="path",
        help="path to the surface calculation directory.",
    )

    parser.add_argument(
        "-epa",
        type=float,
        metavar="epa",
        help="energy per atom of the bulk structure.",
    )

    args = parser.parse_args()

    get_surface_energy(args.epa, args.path)
