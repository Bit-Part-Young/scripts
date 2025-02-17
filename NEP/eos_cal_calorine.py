"""使用 GPUMD & calorine 进行 EOS 计算"""

import argparse

import numpy as np
from ase.io import read, write
from calorine.calculators import CPUNEP, GPUNEP


def eos_cal(
    structure_fn: str,
    model_fn: str,
    num: int = 21,
    output_fn: str = "eos.xyz",
):

    atoms_init = read(structure_fn, format="extxyz")

    # calc = CPUNEP(model_filename=model_fn)
    # 使用 GPUNEP 时，会生成临时文件夹，结束后自动删除，也可选择将其保留
    calc = GPUNEP(model_filename=model_fn)

    for index, volumetric_strain in enumerate(np.linspace(0.9, 1.1, num), start=1):
        atoms = atoms_init.copy()

        atoms.cell *= volumetric_strain

        atoms.calc = calc

        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        stress = atoms.get_stress()

        write(
            output_fn,
            atoms,
            format="extxyz",
            append=True,
        )

    print(f"No. {index} structure cal Done.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="EOS calculation with GPUMD & calorine.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        help="Input structure filename",
    )

    parser.add_argument(
        "model_fn",
        type=str,
        help="NEP model filename",
    )

    parser.add_argument(
        "-n",
        "--num",
        nargs="?",
        const=21,
        default=21,
        type=int,
        help="The number of 0.9 ~ 1.1 volume strain",
    )

    parser.add_argument(
        "-o",
        "--output_fn",
        nargs="?",
        const="eos.xyz",
        default="eos.xyz",
        type=str,
        help="Output structure filename",
    )

    args = parser.parse_args()

    eos_cal(args.structure_fn, args.model_fn, args.num, args.output_fn)
