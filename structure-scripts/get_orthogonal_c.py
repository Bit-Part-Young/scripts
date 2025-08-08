#!/usr/bin/env python3

"""
强制使 c 轴与 a、b 轴垂直（对于生成晶界、界面结构时有用）

reference:
- pymatgen 中的 get_orthogonal_c_slab() 方法
- https://github.com/Ternity/Personal_Scripts/blob/master/ASE/%E5%88%A9%E7%94%A8ase%E5%AE%9E%E7%8E%B0z%E8%BD%B4%E6%96%B9%E5%90%91%E8%BD%AC%E6%AD%A3%E4%BA%A4%E7%9A%84%E6%99%B6%E6%A0%BC%E5%8F%98%E6%8D%A2.md
"""

import argparse

import numpy as np
from ase.io import read, write


def get_orthogonal_c(structure_fn: str = "POSCAR", save_poscar: bool = False):
    """强制使 c 轴与 a、b 轴垂直（对于生成晶界、界面结构时有用）"""

    atoms = read(structure_fn)
    a, b, c = atoms.cell

    # 求 (a, b) 叉积，得垂直于 (a, b) 的向量
    _new_c = np.cross(a, b)
    # 求模(范数), 做归一化,得到单位向量
    _new_c /= np.linalg.norm(_new_c)
    # 求向量 c 在单位向量投影长度, 再 × 单位向量方向
    new_c = np.dot(c, _new_c) * _new_c

    # 根据盒子改变原子坐标
    atoms.set_cell([a, b, new_c], scale_atoms=True)

    if save_poscar:
        output_fn = "orthogonal_c.vasp"
        write(output_fn, atoms, format="vasp", vasp5=True, direct=True)
        print(
            f"\nstructure with orthogonal c (forced) is generated and saved to {output_fn}."
        )
    else:
        print(f"\nstructure with orthogonal c (forced) is generated.")

    return atoms


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate structure with orthogonal c (forced).",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn",
        nargs="?",
        default="POSCAR",
        metavar="structure_fn",
        help="structure filename",
    )
    parser.add_argument("-o", action="store_true", help="write file")

    args = parser.parse_args()

    atoms = get_orthogonal_c(args.structure_fn, args.o)
