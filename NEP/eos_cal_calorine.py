#!/usr/bin/env python3

"""使用 GPUMD & calorine 进行 EOS 计算"""

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

    if os.path.exists(output_fn):
        os.remove(output_fn)

    structure_format = structure_fn.split(".")[-1]
    if structure_format in ["extxyz", "xyz"]:
        atoms_init = read(structure_fn, format="extxyz")
    else:
        atoms_init = read(structure_fn)

    calc = CPUNEP(model_filename=model_fn)
    # 使用 GPUNEP 时，会生成临时文件夹，结束后自动删除，也可选择将其保留
    # calc = GPUNEP(model_filename=model_fn)

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
            "lc": round(atoms.cell.lengths()[0], 5),
            "volume_pa": round(atoms.get_volume() / natoms, 5),
            "energy_pa": round(energy / natoms, 5),
        }

        data_list.append(data_dict)

        write(output_fn, atoms, format="extxyz", append=True)

        print(f"No. {index} structure cal Done.")

    df = pd.DataFrame(data_list)
    df.to_csv("eos_nep.dat", index=False, sep=" ")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="EOS calculation with GPUMD & calorine.",
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
