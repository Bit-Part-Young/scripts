#!/usr/bin/env python3

"""
检查 OUTCAR 中的原子受力收敛情况（使用 ASE 程序）

reference: http://bbs.keinsci.com/thread-19985-1-1.html
"""

import argparse
import os

import numpy as np
from ase.io import read


def grab_outcar_info(outcar_path: str) -> tuple[int, np.ndarray, np.ndarray]:
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
        description="Check the force convergence in OUTCAR.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
    )

    parser.add_argument(
        "-fcc",
        "--ediffg",
        nargs="?",
        const=0.01,
        default=0.01,
        type=float,
        help="force convergence criteria",
    )

    parser.add_argument(
        "outcar_path",
        nargs="?",
        default=".",
        type=str,
        help="OUTCAR path.",
    )
    args = parser.parse_args()

    ediffg = args.ediffg
    outcar_path = args.outcar_path

    natoms, _, force_array = grab_outcar_info(outcar_path=outcar_path)

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
