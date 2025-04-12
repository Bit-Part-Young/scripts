#!/usr/bin/env python3

"""
检查 VASP OUTCAR 文件中的每个离子步原子受力收敛情况（使用 re 正则内置模块）

reference: http://bbs.keinsci.com/thread-19985-1-1.html
"""

import argparse
import os
import re

import numpy as np


def grab_info(outcar_path: str) -> tuple[
    int,
    float,
    np.ndarray,
]:
    """获取原子数、每个离子步中每个原子的受力（向量）"""

    patten_natoms = re.compile(r"\s+\w+\s+NIONS =\s+(\d+)")
    patten_ediffg = re.compile(r"\s+EDIFFG =(\s[-+]?\.\w+[-+]?\d+)")
    patten_force = re.compile(r"\s+total drift:\s+")

    force_list = []
    line_list = []

    outcar_fn = os.path.join(outcar_path, "OUTCAR")
    with open(outcar_fn, "r") as f:
        for index, line in enumerate(f):
            line_list.append(line.strip().split())

            if match := patten_natoms.search(line):
                natoms = int(match.group(1))

            elif match := patten_ediffg.search(line):
                ediffg = float(match.group(1))

            elif patten_force.search(line):
                # 获取一个离子步中的每个原子受力信息
                force_list += line_list[-natoms - 2 : index - 1]

    # 原子位置分量是前 3 列，原子受力分量是后 3 列
    # position_array = np.array(force_list, dtype=float).reshape(-1, natoms, 6)[:, :, :3]
    force_array = np.array(force_list, dtype=float).reshape(-1, natoms, 6)[:, :, 3:]

    return (natoms, ediffg, force_array)


def calcuate_force(
    force_array: np.ndarray,
    force_criteria: float = 0.01,
) -> np.ndarray:
    """计算每个离子步中每个原子的受力（数值），并判断是否达到 EDIFFG 收敛判据"""

    F = np.sqrt(np.power(force_array, 2).sum(axis=-1))

    boolen_array = F > np.abs(force_criteria)

    return boolen_array


def main():
    """主函数"""

    parser = argparse.ArgumentParser(
        description="Check the force convergence of every atom in every ion step in VASP OUTCAR.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
    )

    parser.add_argument(
        "outcar_path",
        nargs="?",
        type=str,
        const=".",
        default=".",
        help="OUTCAR path",
    )

    parser.add_argument(
        "--ediffg",
        nargs="?",
        type=float,
        help="force convergence criteria",
    )

    args = parser.parse_args()

    outcar_path = args.outcar_path
    ediffg = args.ediffg

    if ediffg is None:
        natoms, ediffg, force_array = grab_info(outcar_path=outcar_path)
    else:
        natoms, _, force_array = grab_info(outcar_path=outcar_path)

    ion_steps = force_array.shape[0]
    print(
        f"OUTCAR info: {natoms} atoms, {ion_steps} ion steps, EDIFFG {abs(ediffg)} eV/Å.\n"
    )

    boolen_array = calcuate_force(
        force_array=force_array,
        force_criteria=ediffg,
    )

    # 只输出最后 5 个离子步的受力收敛信息
    print(f"The last 5 ion steps:")
    if ion_steps < 5:
        start = 1
    else:
        boolen_array = boolen_array[-5:, :]
        start = ion_steps - 4

    for step, column in enumerate(boolen_array, start=start):
        # 原子索引从 1 开始
        atom_index = (np.where(column == True)[0] + 1).tolist()

        print(
            f"Step {step}: total {len(atom_index)} atoms force did NOT converge. Index:"
        )

        print(atom_index)


if __name__ == "__main__":
    main()
