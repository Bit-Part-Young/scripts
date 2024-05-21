"""
Check the force convergence in OUTCAR.

reference: http://bbs.keinsci.com/thread-19985-1-1.html

Author: YSL
Version: v0.1
Date: May 21, 2024
"""

import argparse
import re
from pathlib import Path
from typing import Tuple

import numpy as np


def check_outcar(outcar_filepath: Path):
    """Checks whether the OUTCAR file exists."""

    if not outcar_filepath.is_file():
        dashline = "-" * 79
        warning_str = "OUTCAR file does NOT exist! Please check your directory."
        warning_info = "\n".join((dashline, warning_str, dashline))

        raise SystemExit(warning_info)


def grab_info(outcar_filepath: Path) -> Tuple[int, float, np.ndarray, np.ndarray]:
    """Grab the number of atoms, the force convergence criteria,
    position array and force array in all ion steps from the OUTCAR file.

    Returns
    -------
    position_array: numpy 3D array
        The array of atom posistion in all ion steps with unit of angstrom
    force_array: numpy 3D array
        The array of the atom force in all ion steps with unit of eV/Angst
    """

    patten_natoms = re.compile(r"\s+\w+\s+NIONS =\s+(\d+)")
    patten_ediffg = re.compile(r"\s+EDIFFG =(\s[-+]?\.\w+[-+]?\d+)")
    patten_force = re.compile(r"\s+total drift:\s+")

    force_list = []
    line_list = []

    with open(outcar_filepath, "r") as f:
        for index, line in enumerate(f):
            line_list.append(line.strip().split())

            if match := patten_natoms.search(line):
                natoms = int(match.group(1))

            elif match := patten_ediffg.search(line):
                ediffg = float(match.group(1))

            elif patten_force.search(line):
                force_list += line_list[-natoms - 2 : index - 1]

    position_array = np.asfarray(force_list).reshape(-1, natoms, 6)[:, :, :3]
    force_array = np.asfarray(force_list).reshape(-1, natoms, 6)[:, :, 3:]

    return (natoms, ediffg, position_array, force_array)


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
    2) grab the number of atoms, the convergence criteria, and the force
       components of each atom from the OUTCAR file;
    3) calculate the force for each atom;
    4) print the steps and atoms with the force that larger than
       the convergence criteria.
    """

    parser = argparse.ArgumentParser(
        description="Check the force convergence in OUTCAR."
    )
    parser.add_argument(
        "OUTCAR_FILE",
        type=Path,
        help="OUTCAR file",
        default="OUTCAR",
    )
    args = parser.parse_args()

    outcar_filepath = args.OUTCAR_FILE

    check_outcar(outcar_filepath)

    natoms, ediffg, _, force_array = grab_info(outcar_filepath=outcar_filepath)

    print(f"OUTCAR info: {natoms} atoms, {force_array.shape[0]} ion steps.")

    boolen_array = calcuate_force(
        force_array=force_array,
        force_criteria=ediffg,
    )

    for step, column in enumerate(boolen_array, start=1):
        atom_index = (np.where(column == True)[0] + 1).tolist()

        print(
            f"Step {step}: total {len(atom_index)} atoms force did NOT converge, index: {atom_index}."
        )


if __name__ == "__main__":
    main()
