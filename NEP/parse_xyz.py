#!/usr/bin/env python3

"""使用 ASE 模块解析 extxyz 文件并获取 energy, virial/stress, forces 信息"""


import argparse

import numpy as np
from ase.io import read

np.set_printoptions(precision=10, suppress=True)


def parse_xyz(xyz_fn: str, atoms_indices: list[int]):
    """解析 extxyz 文件并获取 energy, virial/stress, forces 信息"""

    atoms_list = read(xyz_fn, index=":", format="extxyz")
    if len(atoms_indices) > len(atoms_list):
        print(f"\nThe atoms indices {atoms_indices} is out of range, Exit!")
        exit(1)

    for i, atoms in enumerate(atoms_list):
        if i not in atoms_indices:
            continue

        try:
            energy = atoms.get_potential_energy()
        except:
            print(f"No energy, forces or stress info in {xyz_fn}, Exit!")
            exit(1)

        print()
        print("-" * 100)

        print(f"The {i} frame info: {atoms.get_chemical_formula()}\n")

        print(f"energy: {round(energy, 5)} eV.\n")

        virial = None
        virial_tag_list = ["virial", "Virial", "virials", "Virials", "stress", "Stress"]
        for virial_tag in virial_tag_list:
            if virial_tag in atoms.info:
                virial = atoms.info.get(virial_tag, "N/A")

                print(f"{virial_tag}:")
                print(virial)
                print()

                break

        if virial is None:
            print("Virial/Stress: N/A.\n")

        forces = None
        forces_tag_list = ["force", "forces"]
        for forces_tag in forces_tag_list:
            if forces_tag in atoms.arrays:
                forces = atoms.arrays.get(forces_tag, "N/A")

                print(f"{forces_tag} (only show first 10 atoms):")
                if forces.shape[0] > 10:
                    print(forces[:10])
                else:
                    print(forces)

                print()

                break

        if forces is None:
            print("Forces: N/A.\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse energy, virial and forces info from extxyz file.",
        epilog="Author: SLY.",
    )

    parser.add_argument("xyz_fn", help="xyz filename")
    parser.add_argument(
        "-ai", "--atoms_indices", nargs="+", type=int, metavar="N", help="atoms indices"
    )

    args = parser.parse_args()

    parse_xyz(args.xyz_fn, args.atoms_indices)
