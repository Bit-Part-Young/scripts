#!/usr/bin/env python3

"""使用 NEP 势函数进行 EOS 计算（cell scaling range: 0.9 ~ 1.1）"""

import argparse
import os

import numpy as np
import pandas as pd
from ase.atoms import Atoms
from ase.io import read, write
from calorine.calculators import CPUNEP, GPUNEP


def eos_cal(
    structure_fn: str, model_fn: str, num: int = 21, output_fn: str = "eos.xyz"
):
    """使用 NEP 势函数进行 EOS 计算（cell scaling range: 0.9 ~ 1.1）"""

    if os.path.exists(output_fn):
        os.remove(output_fn)

    structure_format = structure_fn.split(".")[-1]
    if structure_format in ["extxyz", "xyz"]:
        atoms_init = read(structure_fn, format="extxyz")
    else:
        atoms_init = read(structure_fn)

    # calc = CPUNEP(model_filename=model_fn)
    calc = GPUNEP(model_filename=model_fn)

    data_list = []
    for index, cell_scaling in enumerate(np.linspace(0.9, 1.1, num), start=1):
        atoms: Atoms = atoms_init.copy()

        scaled_positions = atoms_init.get_scaled_positions()
        atoms.cell *= cell_scaling

        atoms.set_scaled_positions(scaled_positions)

        atoms.calc = calc

        energy = atoms.get_potential_energy()

        natoms = len(atoms)
        data_dict = {
            "lc": atoms.cell.lengths()[0],
            "volume_pa": atoms.get_volume() / natoms,
            "energy_pa": energy / natoms,
        }

        data_list.append(data_dict)

        write(output_fn, atoms, format="extxyz", append=True)

        print(f"No. {index} structure cal Done.")

    df = pd.DataFrame(data_list).round(5)
    print(df)

    df.to_csv("eos_nep.dat", index=False, sep=" ")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="EOS calculation with NEP potential for given initialstructure (cell scaling range: 0.9 ~ 1.1).",
        epilog="Author: SLY.",
    )

    parser.add_argument("structure_fn", help="Input structure filename")
    parser.add_argument("model_fn", help="NEP model filename")
    parser.add_argument(
        "-n",
        "--num",
        type=int,
        default=21,
        metavar="N",
        help="The number of 0.9 ~ 1.1 volume strain",
    )
    parser.add_argument(
        "-o",
        "--output_fn",
        type=str,
        default="eos.xyz",
        help="Output structure filename",
    )

    args = parser.parse_args()

    eos_cal(args.structure_fn, args.model_fn, args.num, args.output_fn)
