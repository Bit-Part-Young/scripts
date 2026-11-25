#!/usr/bin/env python3

"""固定（特定）原子层的 x/y/z 方向（VASP selective dynamics）

层编号与 layers_identify.py / layers_count.py 一致（同样的分层算法与容差）。
"""

import argparse

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar

AXES = {"x": 0, "y": 1, "z": 2}


def assign_layers(structure: Structure, axis: str = "z", tol: float = 0.01):
    """按指定笛卡尔轴坐标对原子分层，返回层编号数组（从 1 开始）

    与 layers_identify.py 中的同名函数算法一致：相邻原子坐标差 <= tol 视为同一层；
    若另两个基矢在该轴上分量为 0，则跨越周期边界的原子合并到第 1 层。
    """

    i = AXES[axis]
    eps = 1e-7  # 与 ase get_positions(wrap=True) 一致
    frac = (structure.frac_coords + eps) % 1.0 - eps
    coords = structure.lattice.get_cartesian_coords(frac)[:, i]

    def cluster(x):
        order = np.argsort(x)
        labels = np.empty(len(x), dtype=int)
        labels[order] = np.concatenate(([0], np.cumsum(np.diff(x[order]) > tol))) + 1
        return labels

    labels = cluster(coords)

    cell = structure.lattice.matrix
    others = [j for j in range(3) if j != i]
    if labels.max() > 1 and np.allclose(cell[others, i], 0.0):
        period = cell[i, i]
        if coords.min() + period - coords.max() <= tol:
            coords = np.where(labels == labels.max(), coords - period, coords)
            labels = cluster(coords)

    return labels


def layers_fix(
    structure_fn: str,
    layer_indices: list[int] | None = None,
    fixed_axes: tuple[str, str, str] = ("F", "F", "F"),
    axis: str = "z",
    tol: float = 0.01,
    output_fn: str | None = None,
) -> Structure:
    """固定（特定）原子层的 x/y/z 方向

    layer_indices: 要固定的原子层编号（从 1 开始），None 表示全部原子
    fixed_axes: x/y/z 三个方向的 selective dynamics 标记，F 表示固定，T 表示保持原状
    output_fn: 输出 POSCAR 文件名，None 表示不写文件
    """

    structure = Structure.from_file(structure_fn)
    labels = assign_layers(structure, axis=axis, tol=tol)
    print(f"\nAtoms count: {structure.num_sites}; Layers count: {labels.max()}.\n")

    if layer_indices:
        missing = sorted(set(layer_indices) - set(labels))
        if missing:
            print(f"Warning: layer(s) {missing} not found (valid: 1-{labels.max()}).")
        mask = np.isin(labels, layer_indices)
        print(f"Layer(s) {sorted(set(layer_indices))}: {mask.sum()} atoms fixed.")
    else:
        mask = np.ones(structure.num_sites, dtype=bool)
        print(f"All {mask.sum()} atoms fixed.")

    sd = structure.site_properties.get("selective_dynamics")
    sd = np.ones((structure.num_sites, 3), dtype=bool) if sd is None else np.array(sd, dtype=bool)

    fixed_cols = [k for k, flag in enumerate(fixed_axes) if flag == "F"]
    sd[np.ix_(mask, fixed_cols)] = False
    print(f"Fixed axes (x y z): {' '.join(fixed_axes)}.")

    structure.add_site_property("selective_dynamics", sd.tolist())

    if output_fn:
        Poscar(structure, sort_structure=True).write_file(
            output_fn, direct=True, significant_figures=10
        )
        print(f"\n{output_fn} saved.")

    return structure


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fix atomic layers in x/y/z directions (VASP selective dynamics).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument("structure_fn", nargs="?", default="POSCAR", help="structure filename")
    parser.add_argument(
        "-li",
        "--layer_indices",
        type=int,
        nargs="+",
        metavar="N",
        help="layers to fix (1-based, e.g. 1 2); omit to fix all atoms",
    )
    parser.add_argument(
        "-fa",
        "--fixed_axes",
        nargs=3,
        choices=["F", "T"],
        default=["F", "F", "F"],
        metavar=("X", "Y", "Z"),
        help="selective dynamics flags for x y z: F = fix, T = keep (e.g. F F T)",
    )
    parser.add_argument(
        "--axis", choices=list(AXES), default="z", help="axis perpendicular to atomic layers"
    )
    parser.add_argument(
        "-t", "--tol", type=float, default=0.01, help="tolerance (Å) for atoms in the same layer"
    )
    parser.add_argument(
        "-o",
        "--output",
        nargs="?",
        const="layers_fixed.vasp",
        metavar="FN",
        help="write POSCAR to FN (layers_fixed.vasp if FN omitted)",
    )
    args = parser.parse_args()

    layers_fix(
        args.structure_fn,
        layer_indices=args.layer_indices,
        fixed_axes=args.fixed_axes,
        axis=args.axis,
        tol=args.tol,
        output_fn=args.output,
    )
