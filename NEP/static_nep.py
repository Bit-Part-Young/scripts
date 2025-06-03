#!/usr/bin/env python3

"""NEP 势函数 静态计算"""

import argparse

import numpy as np
from ase.io import read
from ase.units import GPa
from calorine.calculators import CPUNEP


def relax_nep(
    structure_fn: str = "POSCAR",
    potential_fn: str = "nep.txt",
):
    """NEP 势函数 静态计算"""

    atoms = read(structure_fn)
    natoms = len(atoms)
    volume = atoms.get_volume()

    calc = CPUNEP(potential_fn)
    atoms.calc = calc

    energy = atoms.get_potential_energy()
    energy_pa = energy / natoms

    forces = atoms.get_forces()

    stress = atoms.get_stress(voigt=True)
    stress_new = np.array(
        [stress[0], stress[1], stress[2], stress[5], stress[3], stress[4]]
    )

    virial = stress_new * volume

    print(f"Energy (eV/atom): {round(energy_pa, 5)}\n")
    if natoms > 10:
        print(f"Force (eV/Å) (first 10 atoms):")
        print(f"{forces[:10, :]}")
    else:
        print(f"Force (eV/Å):")
        print(f"{forces}")
    print(f"\nStress (GPa): {-stress_new / GPa}")

    print(f"\nVirial (eV/atom): {-virial / natoms}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Static calculation using NEP potential.",
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
