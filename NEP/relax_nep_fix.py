#!/usr/bin/env python3

"""等同于 ISIF=2 以及 F F T的弛豫"""

import os
from ase.io import read, write
import numpy as np

from calorine.calculators import CPUNEP, GPUNEP
from calorine.tools import relax_structure

from ase.atoms import Atoms

from ase.optimize import QuasiNewton, FIRE, LBFGS
from ase.constraints import ExpCellFilter, FixedLine


def get_surface_area(atoms: Atoms):
    a, b, c = atoms.cell.lengths()
    alpha, beta, gamma = atoms.cell.angles()

    surface_area = a * b * np.sin(np.deg2rad(gamma))

    return surface_area


def relax_nep_fix(
    structure_fn: str = "POSCAR",
    potential_fn: str = "nep.txt",
    output_fn: str = "relaxed.xyz",
):

    atoms = read(structure_fn, format="vasp")
    # 删除 constraints，否则可能会报错
    atoms.constraints = None

    surface_area = get_surface_area(atoms)
    coeff = 1.6021766208 * 10**4
    print(f"Surface area: {surface_area} angstrom^2")

    constraint = [FixedLine(atom.index, direction=[0, 0, 1]) for atom in atoms]
    atoms.set_constraint(constraint)

    calc = CPUNEP(potential_fn)
    atoms.calc = calc

    relax_structure(atoms, constant_cell=True)

    energy = atoms.get_potential_energy()
    print(f"Energy: {energy} eV")

    # 删除 constraints，否则会报错
    atoms.constraints = None
    write(output_fn, atoms, format="extxyz", append=True)


if __name__ == "__main__":
    output_fn = "relaxed.xyz"
    if os.path.exists(output_fn):
        os.remove(output_fn)

    for i in range(0, 11):
        print(f"Processing POSCAR.{i}...")
        relax_nep_fix(f"POSCAR.{i}", "nep_TiAlNbMoZrV_20250926.txt", output_fn)
