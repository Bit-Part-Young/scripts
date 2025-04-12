#!/usr/bin/env python3

"""
检查 VASP OUTCAR 文件中的每个离子步原子受力收敛情况（使用 pymatgen 程序）

reference: http://bbs.keinsci.com/thread-19985-1-1.html
"""

import argparse
import os

import numpy as np
from ase.io import read
from pymatgen.io.vasp.outputs import Outcar


def grab_outcar_info(outcar_path: str) -> tuple[int, np.ndarray]:
    """获取原子数、每个离子步中每个原子的受力（向量）"""

    outcar = Outcar(os.path.join(outcar_path, "OUTCAR"))

    # 获取原子位置坐标和受力
    table_data = outcar.read_table_pattern(
        header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
        row_pattern=r"\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        footer_pattern=r"\s--+",
        postprocess=lambda x: float(x),
        last_one_only=False,
    )

    table_array = np.array(table_data)

    # 3D np.ndarray
    forces_array = table_array[:, :, 3:]
    natoms = forces_array.shape[1]

    return (natoms, forces_array)


def calcuate_force(
    force_array: np.ndarray,
    force_criteria: float,
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
        default=".",
        type=str,
        help="OUTCAR path",
    )

    parser.add_argument(
        "--ediffg",
        nargs="?",
        const=0.01,
        default=0.01,
        type=float,
        help="force convergence criteria",
    )

    args = parser.parse_args()

    outcar_path = args.outcar_path
    ediffg = args.ediffg

    natoms, forces_array = grab_outcar_info(outcar_path=outcar_path)

    ion_steps = forces_array.shape[0]

    print(
        f"OUTCAR info: {natoms} atoms, {ion_steps} ion steps, EDIFFG {abs(ediffg)} eV/Å.\n"
    )

    boolen_array = calcuate_force(
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
