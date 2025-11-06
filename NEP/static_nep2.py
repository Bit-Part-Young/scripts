#!/usr/bin/env python3

"""NEP 势函数 静态计算（用于批量计算）"""

import argparse
import os

from ase.io import read, write
from calorine.calculators import CPUNEP, GPUNEP


def static_nep(
    input_xyz_fn: str = "POSCAR",
    potential_fn: str = "nep.txt",
    output_xyz_fn: str = "static_nep.xyz",
):
    """NEP 势函数 静态计算（用于批量计算）"""

    if os.path.exists(output_xyz_fn):
        os.remove(output_xyz_fn)

    atoms_list = read(input_xyz_fn, index=":", format="extxyz")

    # calc = CPUNEP(potential_fn)
    calc = GPUNEP(potential_fn)

    flag = 0
    for atoms in atoms_list:
        # 删除已有的 virial 信息
        if atoms.info.get("virial", None) is not None:
            del atoms.info["virial"]

        atoms.calc = calc

        # 添加改代码，才会计算能量、原子受力和应力；会对相关信息进行覆盖
        atoms.get_potential_energy()

        write(output_xyz_fn, atoms, format="extxyz", append=True)

        flag += 1
        print(f"The {flag} configuration is calculated.")

    print(f"\nBatch static calculation completed. {output_xyz_fn} generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Static calculation using NEP potential (for batch calculations).",
        epilog="Author: SLY.",
    )

    parser.add_argument("input_xyz_fn", help="input xyz filename")
    parser.add_argument("potential_fn", default="nep.txt", help="potential filename")
    parser.add_argument(
        "output_xyz_fn", default="static_nep.xyz", help="output xyz filename"
    )

    args = parser.parse_args()

    static_nep(args.input_xyz_fn, args.potential_fn, args.output_xyz_fn)
