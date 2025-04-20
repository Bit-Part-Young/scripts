#!/usr/bin/env python3

"""NEP 势函数 结构优化"""

import argparse

import numpy as np
from ase.io import read
from calorine.calculators import CPUNEP, GPUNEP
from calorine.tools import relax_structure


def relax_nep(
    structure_fn: str = "POSCAR",
    potential_fn: str = "nep.txt",
):
    """NEP 势函数 结构优化"""

    atoms = read(structure_fn)
    natoms = len(atoms)

    calc = CPUNEP(potential_fn)
    atoms.calc = calc

    energy = atoms.get_potential_energy()
    energy_pa = energy / natoms
    print("Before relaxation:")
    print(f"Lattice constants: {np.round(atoms.cell.lengths(), 5)} angstrom.")
    print(f"Energy_pa: {round(energy_pa, 5)} eV/atom.")

    relax_structure(atoms)

    energy = atoms.get_potential_energy()
    energy_pa = energy / natoms
    print("\nAfter relaxation:")
    print(f"Lattice constants: {np.round(atoms.cell.lengths(), 5)} angstrom.")
    print(f"Energy_pa: {round(energy_pa, 5)} eV/atom.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Relaxation using NEP potential.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        default="POSCAR",
        help="structure filename",
    )

    args = parser.parse_args()

    relax_nep(args.structure_fn)
