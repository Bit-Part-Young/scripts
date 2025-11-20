#!/usr/bin/env python3

"""计算构型/轨迹的净力 net force"""

import argparse

import numpy as np


def parse_xyz(input_xyz_fn: str) -> list:
    """解析 extxyz 文件"""

    with open(input_xyz_fn, "r") as file:
        lines = file.readlines()

    forces_list = []
    i = 0
    while i < len(lines):
        natoms = int(lines[i].strip())
        frame_info = lines[i + 1].strip()

        energy_str = frame_info.split("energy=")[1].split()[0]
        energy = float(energy_str)

        lattice_str = frame_info.split('Lattice="')[1].split('"')[0]
        lattice = np.array(list(map(float, lattice_str.split()))).reshape(3, 3)

        # 只返回原子受力数据
        forces = []
        for j in range(i + 2, i + 2 + natoms):
            atom_info = lines[j].strip().split()
            # 确保有原子受力数据
            if len(atom_info) >= 6:
                fx, fy, fz = (
                    float(atom_info[-3]),
                    float(atom_info[-2]),
                    float(atom_info[-1]),
                )
                forces.append([fx, fy, fz])
            else:
                print(
                    f"\nWarning: some atoms in this extxyz file don't have force info! Exit."
                )
                exit()

        forces_list.append(forces)

        i += 2 + natoms

    return forces_list


def get_net_force(forces: list[list[float]]):
    """计算净力"""

    # 净力: 所有原子的 force 之和（每个原子受力各分量之和，再取模）
    # 与 VASP OUTCAR 中的每个离子步后的 total drift 不同
    net_force = float(np.linalg.norm(np.sum(forces, axis=0)))

    return net_force


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Net force calculation of configuration/trajectory.",
        epilog="Author: SLY.",
    )

    parser.add_argument("xyz_fn", help="extxyz structure/trajectory filename")

    args = parser.parse_args()

    forces_list = parse_xyz(args.xyz_fn)
    net_force_list = [get_net_force(forces) for forces in forces_list]

    print(f"\nNet force value(s):")
    print(net_force_list)
