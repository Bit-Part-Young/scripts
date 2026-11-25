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
    # 兼容 -.1E-01、-0.01、-1E-2 等写法
    patten_ediffg = re.compile(r"\s+EDIFFG\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?)")
    patten_force = re.compile(r"\s+total drift:\s+")

    natoms = None
    ediffg = None
    force_list = []
    line_list = []

    outcar_fn = os.path.join(outcar_path, "OUTCAR")
    with open(outcar_fn, "r") as f:
        for index, line in enumerate(f):
            line_list.append(line.strip().split())

            if match := patten_natoms.search(line):
                natoms = int(match.group(1))

            elif match := patten_ediffg.search(line):
                ediffg = float(match.group(1).replace("D", "E").replace("d", "e"))

            elif patten_force.search(line):
                # 获取一个离子步中的每个原子受力信息
                force_list += line_list[-natoms - 2 : index - 1]

    if natoms is None or not force_list:
        raise SystemExit(f"Error: NIONS or TOTAL-FORCE block not found in {outcar_fn}")

    # 原子位置分量是前 3 列，原子受力分量是后 3 列
    # position_array = np.array(force_list, dtype=float).reshape(-1, natoms, 6)[:, :, :3]
    force_array = np.array(force_list, dtype=float).reshape(-1, natoms, 6)[:, :, 3:]

    return (natoms, ediffg, force_array)


def read_selective_dynamics(poscar_fn: str) -> tuple[np.ndarray, np.ndarray] | None:
    """读取 POSCAR 中的晶格矢量及 Selective dynamics 标记；无 Selective dynamics 时返回 None"""

    if not os.path.isfile(poscar_fn):
        return None

    with open(poscar_fn, "r") as f:
        lines = [line.split() for line in f]

    lattice = np.array([line[:3] for line in lines[2:5]], dtype=float)

    # 第 6 行为元素符号（VASP 5 格式）时，原子数在第 7 行
    index = 5 if re.fullmatch(r"\d+", lines[5][0]) else 6
    natoms = sum(int(n) for n in lines[index])

    if not lines[index + 1][0].lower().startswith("s"):
        return None

    # 跳过 Selective dynamics 行和 Direct/Cartesian 行
    coords = lines[index + 3 : index + 3 + natoms]
    mask = np.array([[flag.upper() == "T" for flag in line[3:6]] for line in coords])

    return (lattice, mask)


def apply_selective_dynamics(
    force_array: np.ndarray,
    lattice: np.ndarray,
    mask: np.ndarray,
) -> np.ndarray:
    """将固定方向上的受力置 0

    Selective dynamics 的 T/F 对应晶格矢量 a、b、c 方向（而非笛卡尔 x、y、z 轴）：
    将受力按晶格矢量分解（F = c @ A），把固定方向的分量置 0 后再转回笛卡尔坐标；
    正交晶胞下等价于直接将对应的笛卡尔分量置 0

    与 VASP 5.4.4 源码 dyna.F 中 SET_SELECTED_FORCES_ZERO 的做法一致（KARDIR/DIRKAR）；
    注意：VASP 使用当前离子步的晶格，这里使用 POSCAR 的晶格，
    仅在 ISIF>=3 且晶胞形状改变时与 VASP 有微小偏差
    """

    components = force_array @ np.linalg.inv(lattice)
    components *= mask

    return components @ lattice


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
        # EDIFFG 缺失、为 0 或为正值（能量判据）时，不能作为力的收敛判据
        if ediffg is None or ediffg >= 0:
            print(
                f"Warning: OUTCAR EDIFFG = {ediffg} is not a force criterion, use 0.01 eV/Å instead (or set --ediffg)."
            )
            ediffg = 0.01
    else:
        natoms, _, force_array = grab_info(outcar_path=outcar_path)

    # 与 VASP 一致：判断收敛时不考虑固定方向上的受力（OUTCAR 中输出的是未处理的原始受力）
    selective = read_selective_dynamics(os.path.join(outcar_path, "POSCAR"))
    if selective is None:
        print("No Selective dynamics in POSCAR, use raw forces.")
    else:
        lattice, mask = selective
        if mask.shape[0] != natoms:
            raise SystemExit("Error: atom number in POSCAR does not match OUTCAR")
        force_array = apply_selective_dynamics(force_array, lattice, mask)
        print(f"Selective dynamics in POSCAR applied ({(~mask).any(axis=1).sum()} atoms with fixed directions).")

    ion_steps = force_array.shape[0]
    print(
        f"OUTCAR info: {natoms} atoms, {ion_steps} ion steps, EDIFFG {abs(ediffg)} eV/Å.\n"
    )

    boolen_array = calculate_force(
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
