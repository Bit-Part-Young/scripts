#!/usr/bin/env python3

"""NEP 势函数结构优化（完全弛豫）"""

import argparse
import os

import numpy as np
from ase.io import read, write
from calorine.calculators import CPUNEP, GPUNEP
from calorine.tools import relax_structure


def relax_nep(structure_fn: str = "POSCAR", potential_fn: str = "nep.txt"):
    """NEP 势函数结构优化（完全弛豫）"""
    static_fn = "static_nep.vasp"
    relaxed_fn = "relaxed_nep.vasp"
    relax_xyz_fn = "relax_nep.xyz"

    if os.path.exists(static_fn):
        os.remove(static_fn)
        os.remove(relaxed_fn)
        os.remove(relax_xyz_fn)

    atoms = read(structure_fn)
    natoms = len(atoms)

    # calc = CPUNEP(potential_fn)
    calc = GPUNEP(potential_fn)

    atoms.calc = calc

    energy = atoms.get_potential_energy()
    energy_pa = energy / natoms
    print("\nBefore relaxation:")
    print(f"Lattice constants: {np.round(atoms.cell.lengths(), 5)} angstrom.")
    print(f"Energy_pa: {round(energy_pa, 5)} eV/atom.")

    # 使用 NEP 静态计算的原子位置与初始导入的构型文件中的会不一致，因此需保存静态计算后的结构
    write(static_fn, atoms, format="vasp", direct=True, sort=True)
    write(relax_xyz_fn, atoms, format="extxyz", append=True)

    relax_structure(atoms)

    energy = atoms.get_potential_energy()
    energy_pa = energy / natoms
    print("\nAfter relaxation:")
    print(f"Lattice constants: {np.round(atoms.cell.lengths(), 5)} angstrom.")
    print(f"Energy_pa: {round(energy_pa, 5)} eV/atom.")

    write(relaxed_fn, atoms, format="vasp", direct=True, sort=True)
    write(relax_xyz_fn, atoms, format="extxyz", append=True)

    print(
        f"\nNEP relaxation completed; Structures saved into {static_fn}, {relaxed_fn} and {relax_xyz_fn}."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Relaxation using NEP potential.")

    parser.add_argument("structure_fn", default="POSCAR", help="structure filename")
    parser.add_argument(
        "potential_fn", default="nep.txt", help="NEP potential filename"
    )

    args = parser.parse_args()

    relax_nep(args.structure_fn, args.potential_fn)
