#!/usr/bin/env python3

"""
生成 dimer 构型
reference: https://github.com/libAtoms/workflow/blob/main/wfl/generate/atoms_and_dimers.py
"""

from ase.atoms import Atoms
from ase.io import write
import numpy as np
import os
import argparse


# [ ] 是否考虑不同元素的 dimer 构型
def dimers_generation(
    element: str,
    bond_length: float = 1.0,
    dimer_nsteps: int = 41,
    dimer_box: str = 15.0,
    dimer_factor_range: tuple = (0.6, 2.5),
    output_fn: str = "dimers.xyz",
):
    """生成 dimer 构型"""

    if os.path.exists(output_fn):
        os.remove(output_fn)

    dimer_distance_min = bond_length * dimer_factor_range[0]
    dimer_distance_max = bond_length * dimer_factor_range[1]
    for dimer_separation in np.linspace(
        dimer_distance_min, dimer_distance_max, dimer_nsteps
    ):
        atoms = Atoms(symbols=f"{element}2", cell=[dimer_box] * 3, pbc=[False] * 3)

        atoms.positions[1, 0] = dimer_separation

        write(output_fn, images=atoms, format="extxyz", append=True)

    print(
        f"Total {dimer_nsteps} dimers configurations from {round(dimer_distance_min, 2)} Å to {round(dimer_distance_max, 2)} Å saved to {output_fn}."
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate dimer configurations.",
        epilog="Author: SLY.",
    )

    parser.add_argument("element", help="element symbol")
    parser.add_argument(
        "bond_length", type=float, nargs="?", default=1.0, help="bond length"
    )
    parser.add_argument(
        "dimer_nsteps", type=int, nargs="?", default=41, help="number of dimer steps"
    )
    parser.add_argument(
        "dimer_box", type=float, nargs="?", default=15.0, help="dimer box"
    )

    args = parser.parse_args()

    dimers_generation(args.element, args.bond_length, args.dimer_nsteps, args.dimer_box)
