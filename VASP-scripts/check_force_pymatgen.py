#!/usr/bin/env python3

"""
检查 VASP OUTCAR 文件中的每个离子步原子受力收敛情况（使用 pymatgen 程序）

reference: http://bbs.keinsci.com/thread-19985-1-1.html
"""

import argparse
import os

import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar


def grab_outcar_info(outcar_path: str) -> tuple[int, float | None, np.ndarray]:
    """获取原子数、每个离子步中每个原子的受力（向量）"""

    outcar = Outcar(os.path.join(outcar_path, "OUTCAR"))

    # 获取 EDIFFG（pymatgen 不解析该参数，使用 read_pattern 读取）
    outcar.read_pattern(
        {"ediffg": r"EDIFFG\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?)"},
        postprocess=lambda x: float(x.replace("D", "E").replace("d", "e")),
        terminate_on_match=True,
    )
    ediffg = outcar.data["ediffg"][0][0] if outcar.data["ediffg"] else None

    # 获取原子位置坐标和受力
    table_data = outcar.read_table_pattern(
        header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
        row_pattern=r"\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        footer_pattern=r"\s--+",
        postprocess=lambda x: float(x),
        last_one_only=False,
    )

    if not table_data:
        raise SystemExit("Error: no TOTAL-FORCE block found in OUTCAR")

    table_array = np.array(table_data)

    # 3D np.ndarray
    forces_array = table_array[:, :, 3:]
    natoms = forces_array.shape[1]

    return (natoms, ediffg, forces_array)


def apply_selective_dynamics(
    outcar_path: str,
    forces_array: np.ndarray,
) -> np.ndarray | None:
    """读取 POSCAR 中的 Selective dynamics，将固定方向上的受力置 0；无 Selective dynamics 时返回 None

    Selective dynamics 的 T/F 对应晶格矢量 a、b、c 方向（而非笛卡尔 x、y、z 轴）：
    将受力按晶格矢量分解，把固定方向的分量置 0 后再转回笛卡尔坐标

    与 VASP 5.4.4 源码 dyna.F 中 SET_SELECTED_FORCES_ZERO 的做法一致（KARDIR/DIRKAR）；
    注意：VASP 使用当前离子步的晶格，这里使用 POSCAR 的晶格，
    仅在 ISIF>=3 且晶胞形状改变时与 VASP 有微小偏差
    """

    poscar_fn = os.path.join(outcar_path, "POSCAR")
    if not os.path.isfile(poscar_fn):
        return None

    poscar = Poscar.from_file(poscar_fn, check_for_potcar=False)
    if poscar.selective_dynamics is None:
        return None

    mask = np.array(poscar.selective_dynamics, dtype=bool)
    if mask.shape[0] != forces_array.shape[1]:
        raise SystemExit("Error: atom number in POSCAR does not match OUTCAR")

    lattice = poscar.structure.lattice.matrix
    components = forces_array @ np.linalg.inv(lattice)
    components *= mask

    return components @ lattice


def calculate_force(
    force_array: np.ndarray,
    force_criteria: float,
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
        default=".",
        type=str,
        help="OUTCAR path",
    )

    parser.add_argument(
        "--ediffg",
        nargs="?",
        const=0.01,
        type=float,
        help="force convergence criteria (default: EDIFFG in OUTCAR)",
    )

    args = parser.parse_args()

    outcar_path = args.outcar_path
    ediffg = args.ediffg

    natoms, outcar_ediffg, forces_array = grab_outcar_info(outcar_path=outcar_path)

    # 与 VASP 一致：默认使用 OUTCAR 中的 EDIFFG；正值为能量判据，不能用于判断受力
    if ediffg is None:
        ediffg = outcar_ediffg
        if ediffg is None or ediffg >= 0:
            print(
                f"Warning: OUTCAR EDIFFG = {ediffg} is not a force criterion, use 0.01 eV/Å instead (or set --ediffg)."
            )
            ediffg = 0.01

    # 与 VASP 一致：判断收敛时不考虑固定方向上的受力（OUTCAR 中输出的是未处理的原始受力）
    constrained = apply_selective_dynamics(outcar_path, forces_array)
    if constrained is None:
        print("No Selective dynamics in POSCAR, use raw forces.")
    else:
        forces_array = constrained
        print("Selective dynamics in POSCAR applied.")

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
