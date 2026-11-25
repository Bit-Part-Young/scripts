#!/usr/bin/env python3

"""统计原子层数目及每层原子数"""

import argparse

import numpy as np
import pandas as pd
from ase.atoms import Atoms
from ase.io import read

AXES = {"x": 0, "y": 1, "z": 2}


def read_structure(structure_fn: str) -> Atoms:
    """读取结构文件（LAMMPS data 需显式指定格式，其余由 ASE 自动识别）"""

    if structure_fn.rsplit(".", 1)[-1] in ("lmp", "data", "lammps-data"):
        return read(structure_fn, format="lammps-data", sort_by_id=False)
    return read(structure_fn)


def assign_layers(atoms: Atoms, axis: str = "z", tol: float = 0.01):
    """按指定笛卡尔轴坐标对原子分层，返回 (层编号数组(从 1 开始), 坐标数组)

    相邻原子坐标差 <= tol 视为同一层；若另两个基矢在该轴上分量为 0（坐标沿该轴周期为 L），
    则跨越周期边界（坐标 ≈0 与 ≈L）的原子合并到第 1 层。
    """

    i = AXES[axis]
    coords = atoms.get_positions(wrap=True)[:, i]

    def cluster(x):
        order = np.argsort(x)
        labels = np.empty(len(x), dtype=int)
        labels[order] = np.concatenate(([0], np.cumsum(np.diff(x[order]) > tol))) + 1
        return labels

    labels = cluster(coords)

    cell = atoms.cell.array
    others = [j for j in range(3) if j != i]
    if atoms.pbc[i] and labels.max() > 1 and np.allclose(cell[others, i], 0.0):
        period = cell[i, i]
        if coords.min() + period - coords.max() <= tol:
            coords = np.where(labels == labels.max(), coords - period, coords)
            labels = cluster(coords)

    return labels, coords


def layers_count(structure_fn: str, axis: str = "z", tol: float = 0.01) -> pd.DataFrame:
    """统计原子层数目及每层原子数"""

    atoms = read_structure(structure_fn)
    labels, coords = assign_layers(atoms, axis=axis, tol=tol)

    df = (
        pd.DataFrame({"layer": labels, axis: coords})
        .groupby("layer")
        .agg(natoms=(axis, "size"), **{f"{axis}(Å)": (axis, "mean")})
        .round(4)
    )

    print(f"Atoms count: {len(atoms)}; Layers count: {len(df)}.\n")
    print(df.to_string())

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count atomic layers and the number of atoms in each layer.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument("structure_fn", nargs="?", default="POSCAR", help="structure filename")
    parser.add_argument(
        "--axis", choices=list(AXES), default="z", help="axis perpendicular to atomic layers"
    )
    parser.add_argument(
        "-t", "--tol", type=float, default=0.01, help="tolerance (Å) for atoms in the same layer"
    )
    args = parser.parse_args()

    layers_count(args.structure_fn, axis=args.axis, tol=args.tol)
