#!/usr/bin/env python3

import argparse

import numpy as np
from ase.atoms import Atom, Atoms
from ase.build import make_supercell
from ase.cell import Cell
from ase.io import write

# HCP 常见类型自间隙原子位置
hcp_sia_dict = {
    "O": [
        [2 / 3, -1 / 3, 0.25],
    ],
    "BO": [
        [2 / 3, -1 / 3, 0.0],
    ],
    "T": [
        [0.0, 0.0, 3 / 8],
    ],
    "BT": [
        [0.0, 0.0, 0.5],
    ],
    "C": [
        [1 / 6, 1 / 6, 0.25],
    ],
    "BC": [
        [0.5, 0.0, 0.0],
    ],
    # S dumbbell 粗略初始位置
    "S": [
        [0.0, 0.0, 0.25],
        [0.0, 0.0, -0.25],
    ],
    # BS dumbbell 粗略初始位置
    "BS": [
        [0.5, 0.0, 0.0],
        [-0.5, 0.0, 0.0],
    ],
}


def hcp_bulk(element: str, a: float, c: float) -> Atoms:
    """HCP 单胞（轴夹角为 60°）"""

    atoms = Atoms(
        symbols=[element] * 2,
        cell=Cell.fromcellpar([a, a, c, 90, 90, 60]),
        scaled_positions=[
            [0.0, 0.0, 0.0],
            [1 / 3, 1 / 3, 0.5],
        ],
        pbc=True,
    )

    return atoms


def hcp_sia_generation(
    element: str,
    a: float,
    c: float,
    sia_type: str,
    size: list[int],
    tranform: bool = False,
    wrap: bool = False,
    output_fn: str = "POSCAR",
):
    """HCP 结构常见类型自间隙原子构型生成"""

    atoms = hcp_bulk(element=element, a=a, c=c)

    supercell = make_supercell(
        prim=atoms,
        P=np.diag(size),
    )

    sia_coords = np.array(hcp_sia_dict[sia_type]) / np.array(size)
    sia_coords_cartesian = sia_coords @ supercell.cell

    for i in range(sia_coords.shape[0]):
        supercell.append(
            Atom(
                symbol=element,
                position=sia_coords_cartesian[i, :],
            )
        )

    # SIA 类型为 S、BS 时，删除坐标为 [0, 0, 0] 的原子
    if sia_type in ["S", "BS"]:
        del supercell[[atom.index for atom in supercell if np.allclose(atom.position, 0.0)]]

    # 将坐标轴夹角 60° 转变成 120°
    if tranform:
        matrix = np.array(
            [
                [1, 0, 0],
                [-1, 1, 0],
                [0, 0, 1],
            ]
        )
        supercell = make_supercell(
            prim=supercell,
            P=matrix,
        )

    if wrap:
        supercell.wrap()

    sia_conc = 1 / (len(supercell) - 1)

    print(f"SIA type: {sia_type}")
    print(f"Supercell size: {size}, natoms: {len(supercell)}")
    print(f"SIA concentration: {sia_conc*100:.8f}%")
    print("\nDefault degree of axes is 60°.")
    print(f"SIA direct coordinates: {sia_coords}")
    print(f"SIA cartesian coordinates: {sia_coords_cartesian}")

    write(output_fn, supercell, format="vasp", direct=True, sort=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate HCP configuration with common SIA type.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument("-e", "--element", type=str, default="Zr", help="Element symbol")

    parser.add_argument("-lc", "--lattice_constant", type=float, nargs=2, help="Lattice constant")

    parser.add_argument("-d", "--dim", type=int, nargs=3, help="x y z dimension")

    parser.add_argument("--sia_type", type=str, help="SIA type")

    parser.add_argument(
        "-t",
        "-tranform",
        action="store_true",
        help="Transform the angle of the axes from 60° to 120°",
    )

    parser.add_argument("-w", "--wrap", action="store_true", help="Wrap the atoms into cell")

    parser.add_argument(
        "-o",
        "--output_fn",
        type=str,
        default="POSCAR",
        help="Output filename",
    )

    args = parser.parse_args()

    hcp_sia_generation(
        element=args.element,
        a=args.lattice_constant[0],
        c=args.lattice_constant[1],
        sia_type=args.sia_type,
        size=args.dim,
        tranform=args.t,
        wrap=args.wrap,
        output_fn=args.output_fn,
    )
