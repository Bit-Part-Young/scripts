#!/usr/bin/env python3

"""
Check the force convergence in OUTCAR with ASE program.

reference: http://bbs.keinsci.com/thread-19985-1-1.html

Author: YSL
Version: v0.1
Date: May 21, 2024
"""

import argparse
import os
import re
from pathlib import Path
from typing import Tuple

import numpy as np
from ase.io import read


def check_outcar(outcar_path: Path):
    """Checks whether the OUTCAR file exists."""

    outcar = os.path.join(outcar_path, "OUTCAR")

    if not outcar.is_file():
        dashline = "-" * 79
        warning_str = (
            "OUTCAR file does NOT exist! Please check your directory."
        )
        warning_info = "\n".join((dashline, warning_str, dashline))

        raise SystemExit(warning_info)


def grab_info(outcar_path: Path) -> Tuple[int, np.ndarray, np.ndarray]:
    """Grab the number of atoms, the force convergence criteria,
    position array and force array in all ion steps from the OUTCAR file.

    Returns
    -------
    position_array: numpy 3D array
        The array of atom posistion in all ion steps with unit of angstrom
    force_array: numpy 3D array
        The array of the atom force in all ion steps with unit of eV/Angst
    """

    outcar = os.path.join(outcar_path, "OUTCAR")
    atoms_list = read(outcar, index=":")
    natoms = len(atoms_list[0])

    position_array = np.array([atoms.get_positions() for atoms in atoms_list])
    force_array = np.array([atoms.get_forces() for atoms in atoms_list])

    return (natoms, position_array, force_array)


def calcuate_force(
    force_array: np.ndarray,
    force_criteria: float,
):
    """Calculate the force for each atom and compare with the force convergence criteria.

    Parameters
    ----------
    force_array: numpy 3D array
        The array of atom posistion in all ion steps with unit of angstrom

    Returns
    -------
    Boolean array
        The Boolean array that True represents the force of the atom
    does NOT converged.
    """

    F = np.sqrt(np.power(force_array, 2).sum(axis=-1))

    boolen_array = F > np.abs(force_criteria)

    return boolen_array


def main():
    """Workflow
    1) check whether the OUTCAR file exists or not;
    2) grab the number of atoms and the force components
       of each atom from the OUTCAR file;
    3) calculate the force for each atom;
    4) print the steps and atoms with the force
       that larger than the convergence criteria.
    """

    parser = argparse.ArgumentParser(
        description="Check the force convergence in OUTCAR."
    )
    parser.add_argument(
        "-ediffg",
        type=float,
        help="force convergence criteria. Default: 0.01.",
        default=0.01,
    )
    parser.add_argument(
        "OUTCAR_PATH",
        type=Path,
        help="OUTCAR file path. Default: .",
        default="OUTCAR",
    )
    args = parser.parse_args()

    outcar_path = args.OUTCAR_PATH
    ediffg = args.ediffg

    check_outcar(outcar_path)

    natoms, _, force_array = grab_info(outcar_path=outcar_path)

    print(f"OUTCAR info: {natoms} atoms, {force_array.shape[0]} ion steps.\n")

    boolen_array = calcuate_force(
        force_array=force_array,
        force_criteria=float(ediffg),
    )

    for step, column in enumerate(boolen_array, start=1):
        atom_index = (np.where(column == True)[0] + 1).tolist()

        print(
            f"Step {step}: Total {len(atom_index)} atoms force did NOT converge, Index:"
        )
        print(atom_index)


if __name__ == "__main__":
    main()
