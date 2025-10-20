#!/usr/bin/env python3

"""给定 EOS 构型文件，使用 NEP 势函数进行 EOS 计算"""

import argparse

import numpy as np
import pandas as pd
from ase.io import read
from calorine.calculators import CPUNEP, GPUNEP


def eos_nep(structure_fn: str = "eos.xyz", potential_fn: str = "nep.txt"):
    """给定 EOS 构型文件，使用 NEP 势函数进行 EOS 计算"""

    atoms_list = read(structure_fn, index=":", format="extxyz")

    # calc = CPUNEP(model_filename=potential_fn)
    calc = GPUNEP(model_filename=potential_fn)

    data_list = []
    for atoms in atoms_list:
        atoms.calc = calc

        natoms = len(atoms)
        volume = atoms.get_volume()
        energy = atoms.get_potential_energy()

        cellsize = int(np.ceil((natoms / 2) ** (1 / 3)))
        lc = atoms.cell.lengths()[0] / cellsize
        volume_pa = volume / natoms
        energy_pa = energy / natoms

        data_dict = {"lc": lc, "volume_pa": volume_pa, "energy_pa": energy_pa}
        data_list.append(data_dict)

    df = pd.DataFrame(data_list).round(5)
    print(df)

    df.to_csv("eos_nep.dat", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="EOS calculation with NEP potential for given EOS xyz file.",
        epilog="Author: SLY.",
    )

    parser.add_argument("structure_fn", default="eos.xyz", help="eos xyz filename")
    parser.add_argument("potential_fn", default="nep.txt", help="nep model filename")

    args = parser.parse_args()

    eos_nep(args.structure_fn, args.potential_fn)
