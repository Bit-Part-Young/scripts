#!/usr/bin/env python3

"""
检查 VASP OUTCAR 文件中的每个离子步原子受力收敛情况（使用 ASE 程序）

reference: http://bbs.keinsci.com/thread-19985-1-1.html
"""

# [ ] 添加力未收敛中的最大值？

import argparse
import os
import re

import numpy as np
from ase.io import read


def grab_outcar_info(outcar_path: str) -> tuple[int, np.ndarray, bool]:
    """获取原子数、每个离子步中每个原子的受力（向量）"""

    outcar = os.path.join(outcar_path, "OUTCAR")
    atoms_list = read(outcar, index=":", format="vasp-out")
    natoms = len(atoms_list[0])

    # ASE 会读取 OUTCAR 同目录下的 CONTCAR（非空时优先）或 POSCAR 中的 Selective dynamics，
    # 转为 FixAtoms / FixScaled 约束；get_forces() 默认会将原子固定轴对应的受力直接设置为 0，
    # 与 VASP 判断收敛时的做法一致。若两个文件都不存在，则为原始受力
    # 注意：VASP（dyna.F 中 SET_SELECTED_FORCES_ZERO）用倒格矢分解受力 c_j = F·b_j，
    # 置 0 后再 F' = sum_j c_j·a_j；ASE 的 FixScaled 则用 a_j·F 分解，
    # 因此非正交晶胞中部分方向固定（如 T T F）的原子，结果与 VASP 不同；
    # 正交晶胞、全部固定（F F F）或全部放开的原子则完全一致
    # 3D np.ndarray
    forces_array = np.array([atoms.get_forces() for atoms in atoms_list])
    constrained = bool(atoms_list[0].constraints)

    return (natoms, forces_array, constrained)


def grab_ediffg(outcar_path: str) -> float | None:
    """从 OUTCAR 中读取 EDIFFG（ASE 不解析该参数）"""

    pattern = re.compile(r"EDIFFG\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?)")
    with open(os.path.join(outcar_path, "OUTCAR"), "r") as f:
        for line in f:
            if match := pattern.search(line):
                return float(match.group(1).replace("D", "E").replace("d", "e"))

    return None


def calculate_force(
    force_array: np.ndarray,
    force_criteria: float = 0.01,
) -> np.ndarray:
    """计算每个离子步中每个原子的受力（数值），并判断是否达到 EDIFFG 收敛判据"""

    F = np.sqrt(np.power(force_array, 2).sum(axis=-1))

    # 与 VASP 源码 main.F 一致：受力模长 > |EDIFFG| 视为未收敛
    # 注意：ISIF>=3 时 VASP 还要求应力满足 EDIFFG，这里只检查受力
    boolen_array = F > np.abs(force_criteria)

    return boolen_array


def main():
    """主函数"""

    parser = argparse.ArgumentParser(
        description="Check the force convergence of every atom in every ion step in VASP OUTCAR.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
        epilog="Author: SLY.",
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
        const=0.01,
        help="force convergence criteria (default: EDIFFG in OUTCAR)",
    )

    args = parser.parse_args()

    outcar_path = args.outcar_path
    ediffg = args.ediffg

    natoms, forces_array, constrained = grab_outcar_info(outcar_path=outcar_path)

    # 与 VASP 一致：默认使用 OUTCAR 中的 EDIFFG；正值为能量判据，不能用于判断受力
    if ediffg is None:
        ediffg = grab_ediffg(outcar_path)
        if ediffg is None or ediffg >= 0:
            print(
                f"Warning: OUTCAR EDIFFG = {ediffg} is not a force criterion, use 0.01 eV/Å instead (or set --ediffg)."
            )
            ediffg = 0.01

    if constrained:
        print("Selective dynamics in CONTCAR/POSCAR applied by ASE.")
    else:
        print("No Selective dynamics found by ASE, use raw forces.")

    ion_steps = forces_array.shape[0]

    print(
        f"OUTCAR info: {natoms} atoms, {ion_steps} ion steps, EDIFFG {abs(ediffg)} eV/Å.\n"
    )

    boolen_array = calculate_force(
        force_array=forces_array,
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
            f"Step {step}: Total {len(atom_index)} atoms force did NOT converge. Index:"
        )
        print(atom_index)


if __name__ == "__main__":
    main()
