#!/usr/bin/env python3

"""在指定轴（默认 z 轴）顶部、底部、顶部底部两端添加真空层（支持非正交胞，保持晶格角度不变）

真空沿该轴对应的晶格基矢方向添加，使该基矢在该笛卡尔轴上的分量增加 vacuum（两端各 vacuum）。
注意：原子应位于晶胞内且不跨越该方向的周期边界，否则真空层会把原子层切开。
"""

import argparse

import numpy as np
from ase.atoms import Atoms
from ase.io import read, write

AXES = {"x": 0, "y": 1, "z": 2}


def add_vacuum(
    structure_fn: str,
    vacuum: float = 10.0,
    mode: str = "both",
    axis: str = "z",
    output_fn: str | None = None,
) -> Atoms:
    """在指定轴顶部（top）、底部（bottom）或两端（both）添加厚度为 vacuum (Å) 的真空层"""

    if mode not in ("top", "bottom", "both"):
        raise ValueError(f"Invalid mode: {mode}. Must be 'top', 'bottom', or 'both'.")

    atoms = read(structure_fn)
    i = AXES[axis]
    cell = atoms.cell.array.copy()

    if abs(cell[i, i]) < 1e-8:
        raise ValueError(f"Lattice vector {i + 1} has no component along {axis}.")

    scaled = atoms.get_scaled_positions(wrap=False)[:, i]
    if scaled.min() < -1e-6 or scaled.max() >= 1.0 - 1e-6:
        print(f"Warning: some atoms lie outside the cell along {axis}; check the result.")

    # 沿原基矢方向、且在该轴上分量恰为 vacuum 的向量（保持晶格角度不变）
    vacuum_vector = cell[i] * vacuum / cell[i, i]

    cell[i] += vacuum_vector * (2 if mode == "both" else 1)
    atoms.set_cell(cell, scale_atoms=False)
    if mode in ("bottom", "both"):
        atoms.translate(vacuum_vector)

    print(f"\nAdded {vacuum} Å vacuum along {axis} ({mode}).")

    if output_fn:
        write(output_fn, atoms, format="vasp", direct=True)
        print(f"{output_fn} saved.")

    return atoms


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add vacuum to the top, bottom, or both ends along the specified axis "
        "(lattice angles are preserved).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument("structure_fn", nargs="?", default="POSCAR", help="structure filename")
    parser.add_argument(
        "vacuum", nargs="?", type=float, default=10.0, help="vacuum thickness (Å) on each side"
    )
    parser.add_argument(
        "mode", nargs="?", choices=["top", "bottom", "both"], default="both", help="where to add vacuum"
    )
    parser.add_argument("--axis", choices=list(AXES), default="z", help="axis to add vacuum along")
    parser.add_argument(
        "-o",
        "--output",
        nargs="?",
        const="vacuum.vasp",
        metavar="FN",
        help="write POSCAR to FN (vacuum.vasp if FN omitted)",
    )
    args = parser.parse_args()

    add_vacuum(args.structure_fn, args.vacuum, args.mode, axis=args.axis, output_fn=args.output)
