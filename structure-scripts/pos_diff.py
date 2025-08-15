#!/usr/bin/env python3

"""比较两个构型之间的原子坐标变化/差异（主要为弛豫/变形前后）"""

import argparse
import pprint

import numpy as np
import pandas as pd
from ase.io import read


def pos_diff(
    structure1_fn: str,
    structure2_fn: str,
    wrap: bool = False,
    atom_indices: int | list[int] | None = None,
) -> np.ndarray:
    """比较两个构型之间的原子坐标变化/差异（主要为弛豫/变形前后）"""

    atoms1 = read(structure1_fn)
    atoms2 = read(structure2_fn)

    structure1_info = {
        "natom": len(atoms1),
        "cell": atoms1.cell[:].round(5),
        "volume": float(round(atoms1.get_volume(), 1)),
    }

    structure2_info = {
        "natom": len(atoms2),
        "cell": atoms2.cell[:].round(5),
        "volume": float(round(atoms2.get_volume(), 1)),
    }

    print(f"{structure1_fn} info:")
    pprint.pprint(structure1_info)
    print(f"\n{structure2_fn} info:")
    pprint.pprint(structure2_info)

    vol_diff = (
        (structure1_info["volume"] - structure2_info["volume"])
        / structure2_info["volume"]
    ) * 100

    print(f"\nVolume difference: {vol_diff:.5f}%")

    print(f"\nCoordinate difference between {structure1_fn} with {structure2_fn}:\n")

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

    if atom_indices is None:
        df_selected = df.copy()
    elif isinstance(atom_indices, int) or isinstance(atom_indices, list):
        df_selected = df.iloc[atom_indices]

    print(df_selected)

    # 输出每列的绝对值最大值，转置
    print("\nMax absolute change value:")
    print(df_selected.abs().max().to_frame().T)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Compare the cell & positions difference between two structures (mainly for relaxation/deformation).",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
    )

    parser.add_argument(
        "structure1_fn",
        nargs="?",
        type=str,
        default="CONTCAR",
        metavar="structure1_fn",
        help="structure1 filename after relaxation/deformation",
    )

    parser.add_argument(
        "structure2_fn",
        nargs="?",
        type=str,
        default="POSCAR",
        metavar="structure2_fn",
        help="structure2 filename before relaxation/deformation",
    )

    parser.add_argument(
        "-w", "--wrap", action="store_true", help="wrap atoms into cell"
    )

    parser.add_argument(
        "-ai",
        "--atom_indices",
        nargs="*",
        type=int,
        metavar="atom_indices",
        help="atom indices",
    )

    args = parser.parse_args()

    structure1_fn = args.structure1_fn
    structure2_fn = args.structure2_fn
    atom_indices = args.atom_indices
    wrap = args.wrap

    pos_diff(
        structure1_fn=structure1_fn,
        structure2_fn=structure2_fn,
        atom_indices=atom_indices,
        wrap=wrap,
    )
