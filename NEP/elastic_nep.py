#!/usr/bin/env python3

"""NEP 势函数 弹性常数计算"""

import argparse

import numpy as np
from ase.io import read
from calorine.calculators import CPUNEP, GPUNEP
from calorine.tools import get_elastic_stiffness_tensor, relax_structure

np.set_printoptions(precision=1, suppress=True)


def elastic_nep(
    structure_fn: str = "POSCAR",
    potential_fn: str = "nep.txt",
):
    """NEP 势函数 弹性常数计算"""

    atoms = read(structure_fn)

    calc = CPUNEP(potential_fn)

    atoms.calc = calc

    # 结构优化
    relax_structure(atoms)

    # 弹性常数计算
    cij = get_elastic_stiffness_tensor(atoms)

    print(cij)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate elastic constants using NEP potential.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        help="structure filename",
        default="POSCAR",
    )

    args = parser.parse_args()

    elastic_nep(args.structure_fn)
