#!/usr/bin/env python3

"""查看两个构型之间的原子坐标变化/差异（主要为弛豫前后）"""

import argparse

import numpy as np
import pandas as pd
from ase.io import read


def pos_diff(
    structure1_fn: str,
    structure2_fn: str,
    wrap: bool = False,
    atom_index: int | list[int] | None = None,
) -> np.ndarray:
    """查看两个构型之间的原子坐标变化/差异（主要为弛豫前后）"""

    atoms1 = read(structure1_fn)
    atoms2 = read(structure2_fn)

    structure1_info = {
        "natom": len(atoms1),
        "cell": atoms1.cell,
        "volume": round(atoms1.get_volume(), 1),
    }

    structure2_info = {
        "natom": len(atoms2),
        "cell": atoms2.cell,
        "volume": round(atoms2.get_volume(), 1),
    }

    print(f"{structure1_fn} info: {structure1_info}")
    print(f"{structure2_fn} info: {structure2_info}")

    print(
        f"\nCoordinate difference between {structure1_fn} with {structure2_fn}:\n"
    )

    position1 = atoms1.get_positions(wrap=wrap)
    position2 = atoms2.get_positions(wrap=wrap)

    scaled_positions1 = atoms1.get_scaled_positions(wrap=wrap)
    scaled_positions2 = atoms2.get_scaled_positions(wrap=wrap)

    diff = position1 - position2
    scaled_diff = scaled_positions1 - scaled_positions2

    data = {
        "x": diff[:, 0],
        "y": diff[:, 1],
        "z": diff[:, 2],
        "xs": scaled_diff[:, 0],
        "ys": scaled_diff[:, 1],
        "zs": scaled_diff[:, 2],
    }

    df = pd.DataFrame(data).round(5)

    if atom_index is None:
        print(df)
    elif isinstance(atom_index, int) or isinstance(atom_index, list):
        print(df.iloc[atom_index])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Show atomic positions difference between two structures (for relaxation).",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
    )

    parser.add_argument(
        "structure1_fn",
        nargs="?",
        type=str,
        default="CONTCAR",
        help="Structure1 filename (After relaxation).",
    )

    parser.add_argument(
        "structure2_fn",
        nargs="?",
        type=str,
        default="POSCAR",
        help="Structure2 filename (Before relaxation).",
    )

    parser.add_argument(
        "-w",
        "--wrap",
        type=int,
        nargs="?",
        const=1,
        choices=[0, 1],
        default=None,
        help="Whether to wrap atomic positions into the cell.",
    )

    parser.add_argument(
        "-ai",
        "--atom_index",
        nargs="*",
        type=int,
        help="Atom index.",
    )

    args = parser.parse_args()

    structure1_fn = args.structure1_fn
    structure2_fn = args.structure2_fn
    atom_index = args.atom_index
    wrap = args.wrap

    pos_diff(
        structure1_fn=structure1_fn,
        structure2_fn=structure2_fn,
        atom_index=atom_index,
        wrap=wrap,
    )
