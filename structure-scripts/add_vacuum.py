#!/usr/bin/env python3

"""在指定轴（默认 z 轴）顶部、底部、顶部底部两端添加真空层（支持非正交胞，保持晶格角度不变）"""

import argparse

import numpy as np
from ase.atoms import Atoms
from ase.io import read, write


def add_vacuum(
    structure_fn: str,
    vacuum: float,
    mode: str,
    axis: str = "z",
    save_poscar: bool = False,
) -> Atoms:
    """在指定轴（默认 z 轴）顶部、底部、顶部底部两端添加真空层（支持非正交胞，保持晶格角度不变）"""

    atoms = read(structure_fn, format="vasp")

    cell = atoms.get_cell()

    if axis == "z":
        axis_index = 2
    elif axis == "y":
        axis_index = 1
    elif axis == "x":
        axis_index = 0
    else:
        raise ValueError(f"Invalid axis: {axis}. Must be 'x', 'y', or 'z'.")

    # 指定轴的晶格基矢
    lattice_vector = cell[axis_index]
    # 归一化得到单位基矢
    unit_vector = lattice_vector / np.linalg.norm(lattice_vector)

    # 添加的真空层对应的向量，以保持晶格角度不变
    vacuum_vector = vacuum * unit_vector

    # 确保指定轴的分量增量值 为 vacuum
    projection = vacuum_vector[axis_index]
    if abs(projection - vacuum) > 1e-6:
        scale_factor = vacuum / projection
        vacuum_vector *= scale_factor

    # 在顶部添加真空层
    if mode == "top":
        cell[axis_index] += vacuum_vector
    # 在底部添加真空层
    elif mode == "bottom":
        cell[axis_index] += vacuum_vector
        positions = atoms.get_positions()
        positions += vacuum_vector
        atoms.set_positions(positions)
    # 在两端添加真空层
    elif mode == "both":
        cell[axis_index] += 2 * vacuum_vector
        positions = atoms.get_positions()
        positions += vacuum_vector
        atoms.set_positions(positions)
    else:
        raise ValueError(f"Invalid mode: {mode}. Must be 'top', 'bottom', or 'both'.")

    atoms.set_cell(cell)

    if save_poscar:
        output_fn = "vacuum.vasp"
        write(output_fn, atoms, format="vasp", vasp5=True, direct=True)
        print(f"\nAdded {vacuum} Å vacuum to {mode} side and saved to {output_fn}.")
    else:
        print(f"\nAdded {vacuum} Å vacuum to {mode} side.")

    return atoms


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Add vacuum to the top, bottom, or both ends in specified axis.",
        epilog="Author: SLY.",
    )
    parser.add_argument(
        "structure_fn", nargs="?", default="POSCAR", help="structure filename"
    )
    parser.add_argument(
        "vacuum",
        type=float,
        nargs="?",
        default=10.0,
        help="vacuum thickness (default: 10.0)",
    )
    parser.add_argument(
        "mode",
        choices=["top", "bottom", "both"],
        nargs="?",
        default="both",
        help="mode to add vacuum (default: both)",
    )

    parser.add_argument(
        "--axis",
        choices=["x", "y", "z"],
        default="z",
        metavar="STR",
        help="axis to add vacuum (default: z)",
    )

    parser.add_argument("-o", action="store_true", help="write file")

    args = parser.parse_args()

    add_vacuum(args.structure_fn, args.vacuum, args.mode, args.axis, args.o)
