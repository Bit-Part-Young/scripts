#!/usr/bin/env python3

"""NEP 势函数结构优化（完全弛豫；可单个或批量计算）"""

import argparse
import os

import numpy as np
from ase.atoms import Atoms
from ase.io import read, write
from calorine.calculators import CPUNEP, GPUNEP
from calorine.tools import relax_structure


def relax_nep_single(atoms: Atoms, potential_fn: str = "nep.txt", output: bool = False):

    # calc = CPUNEP(potential_fn)
    calc = GPUNEP(potential_fn)
    atoms.calc = calc

    natoms = len(atoms)
    energy = atoms.get_potential_energy()
    energy_pa = energy / natoms
    forces = atoms.get_forces()

    # copy() 方法只会拷贝 Atoms，其 calculator 不会被拷贝
    atoms_static = atoms.copy()
    atoms_static.info["energy"] = energy
    atoms_static.arrays["forces"] = forces

    if output:
        print("\nBefore relaxation:")
        print(f"Lattice constants: {np.round(atoms.cell.lengths(), 5)} angstrom.")
        print(f"Energy_pa: {round(energy_pa, 5)} eV/atom.")

    relax_structure(atoms, minimizer="fire")

    energy = atoms.get_potential_energy()
    energy_pa = energy / natoms
    forces = atoms.get_forces()

    atoms_relaxed = atoms.copy()
    atoms_relaxed.info["energy"] = energy
    atoms_relaxed.arrays["forces"] = forces

    if output:
        print("\nAfter relaxation:")
        print(f"Lattice constants: {np.round(atoms.cell.lengths(), 5)} angstrom.")
        print(f"Energy_pa: {round(energy_pa, 5)} eV/atom.")

    return atoms_static, atoms_relaxed


def relax_nep(structure_fn: str = "POSCAR", potential_fn: str = "nep.txt"):
    """NEP 势函数结构优化（完全弛豫）"""

    input_fn_basename = os.path.basename(structure_fn)
    input_format = input_fn_basename.split(".")[-1]

    static_fn = "static_nep.vasp"
    relaxed_fn = "relaxed_nep.vasp"
    relax_xyz_fn = "relax_nep.xyz"

    if os.path.exists(static_fn):
        os.remove(static_fn)
        os.remove(relaxed_fn)
        os.remove(relax_xyz_fn)

    if input_format in ["xyz", "extxyz"]:
        atoms_list: list[Atoms] = read(structure_fn, index=":", format="extxyz")

        static_fn = "static_nep.xyz"
        relaxed_fn = "relaxed_nep.xyz"
        if os.path.exists(static_fn):
            os.remove(static_fn)
            os.remove(relaxed_fn)

        for atoms in atoms_list:
            atoms_static, atoms_relaxed = relax_nep_single(
                atoms, potential_fn, output=False
            )
            write(static_fn, atoms_static, format="extxyz", append=True)
            write(relaxed_fn, atoms_relaxed, format="extxyz", append=True)
    else:
        atoms = read(structure_fn)

        static_fn = "static_nep.vasp"
        relaxed_fn = "relaxed_nep.vasp"
        relax_xyz_fn = "relax_nep.xyz"
        if os.path.exists(static_fn):
            os.remove(static_fn)
            os.remove(relaxed_fn)
            os.remove(relax_xyz_fn)

        atoms_static, atoms_relaxed = relax_nep_single(atoms, potential_fn, output=True)

        # 使用 NEP 静态计算的原子位置与初始导入的构型文件中的会不一致，因此需保存静态计算后的结构
        write(static_fn, atoms_static, format="vasp", direct=True, sort=True)
        write(relaxed_fn, atoms_relaxed, format="vasp", direct=True, sort=True)
        # 讲静态、弛豫后的构型数据一起保存到 xyz 文件中
        write(relax_xyz_fn, [atoms_static, atoms_relaxed], format="extxyz", append=True)

    print(
        f"\nNEP relaxation completed; Structures saved into {static_fn}, {relaxed_fn} and {relax_xyz_fn}."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Relaxation calculation using NEP potential (POSCAR for single calculation & check; xyz for batch calculation).",
        epilog="Author: SLY.",
    )

    parser.add_argument("structure_fn", default="POSCAR", help="structure filename")
    parser.add_argument(
        "potential_fn", default="nep.txt", help="NEP potential filename"
    )

    args = parser.parse_args()

    relax_nep(args.structure_fn, args.potential_fn)
