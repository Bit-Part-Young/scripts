#!/usr/bin/env python3

"""生成用于 BCC/HCP 结构纯元素 VASP 空位迁移能计算的含空位的初始和末态结构"""

import argparse
import os

import numpy as np
from ase.build import bulk
from ase.geometry.geometry import get_distances
from ase.io import write


def vacancy_migration_generation(symbol: str, a: float, size: int = 3):

    bulk_atoms = bulk(symbol, a=a, cubic=True)
    bulk_atoms *= size

    # BCC 结构，空位原子分数位置刚好是 [0.5, 0.5, 0.5]；FCC 结构不是
    # 将空位放置在 cell 中心附近
    D, D_len = get_distances(
        np.diag(bulk_atoms.cell) / 2,
        bulk_atoms.positions,
        bulk_atoms.cell,
        bulk_atoms.pbc,
    )

    vacancy_index = D_len.argmin()
    vacancy_position = bulk_atoms.positions[vacancy_index]
    vacancy_atoms_initial = bulk_atoms.copy()
    del vacancy_atoms_initial[vacancy_index]

    # BCC 结构，迁移方向 <111>；FCC 结构，迁移方向 <110>
    # 识别空位两个相反的最近邻原子索引 identify two opposing nearest neighbours of the vacancy
    D, D_len = get_distances(
        vacancy_position,
        vacancy_atoms_initial.positions,
        vacancy_atoms_initial.cell,
        vacancy_atoms_initial.pbc,
    )
    D = D[0, :]
    D_len = D_len[0, :]

    nn_mask = np.abs(D_len - D_len.min()) < 1e-8
    # 统计空位原子位置的最近邻原子个数，以判断迁移方向: BCC 8 <111>; FCC 12 <110>
    # print(len(nn_mask.nonzero()[0]))
    i1 = nn_mask.nonzero()[0][0]
    # 空位指向各原子的向量加上 D[i1]，最小化即为与 i1 相反方向的最近邻
    i2 = ((D + D[i1]) ** 2).sum(axis=1).argmin()

    # 将近邻原子放置在空位位置，得到末态结构
    vacancy_atoms_final = vacancy_atoms_initial.copy()
    vacancy_atoms_final.positions[i1] = vacancy_position

    write("initial.vasp", vacancy_atoms_initial, format="vasp", direct=True)
    write("final.vasp", vacancy_atoms_final, format="vasp", direct=True)

    xyz_fn = f"vacancy_migration_{symbol}.xyz"
    if os.path.exists(xyz_fn):
        os.remove(xyz_fn)

    write(
        xyz_fn,
        [vacancy_atoms_initial, vacancy_atoms_final],
        format="extxyz",
        append=True,
    )

    print(f"\ninitial.vasp, final.vasp, and {xyz_fn} generated.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate initial and final configurations for VASP vacancy migration energy calculation (currently only support BCC/HCP structure).",
        epilog="Author: SLY.",
    )

    parser.add_argument("symbol", help="element symbol. eg. Nb")
    parser.add_argument("a", type=float, help="lattice constant. eg. 3.30")
    parser.add_argument("size", type=int, default=3, help="cell size. eg. 3")

    args = parser.parse_args()

    vacancy_migration_generation(args.symbol, args.a, args.size)
