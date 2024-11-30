#!/usr/bin/env python3

"""查看两个构型之间的变化/差异（主要为弛豫前后）"""

import sys

import numpy as np
import pandas as pd
from ase.io import read


def pos_diff(
    structure1_fn: str,
    structure2_fn: str,
    atom_sequence: int | list[int] | None = None,
) -> np.ndarray:
    """查看两个构型之间的变化/差异（主要为弛豫前后）"""

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

    position1 = atoms1.get_positions(wrap=True)
    position2 = atoms2.get_positions(wrap=True)

    scaled_positions1 = atoms1.get_scaled_positions(wrap=True)
    scaled_positions2 = atoms2.get_scaled_positions(wrap=True)

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

    if atom_sequence is None:
        print(df)
    elif isinstance(atom_sequence, int) or isinstance(atom_sequence, list):
        print(df.iloc[atom_sequence])


if __name__ == "__main__":
    structure1_fn = sys.argv[1]
    structure2_fn = sys.argv[2]

    atoms1 = read(structure1_fn)
    atoms2 = read(structure2_fn)

    if len(sys.argv) == 3:
        pos_diff(
            structure1_fn=structure1_fn,
            structure2_fn=structure2_fn,
        )
    elif len(sys.argv) > 3:
        atom_sequence = [int(arg) for arg in sys.argv[3:]]

        pos_diff(
            structure1_fn=structure1_fn,
            structure2_fn=structure2_fn,
            atom_sequence=atom_sequence,
        )
