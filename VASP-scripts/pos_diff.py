#!/usr/bin/env python

"""查看两个构型之间的变化/差异（主要为弛豫前后）"""

import sys

import numpy as np
from ase.atoms import Atoms
from ase.io import read


def pos_diff(
    atoms1: Atoms,
    atoms2: Atoms,
    num: int = None,
) -> np.ndarray:
    """查看两个构型之间的变化/差异（主要为弛豫前后）"""

    position1 = atoms1.get_positions(wrap=True)
    position2 = atoms2.get_positions(wrap=True)

    diff: np.ndarray
    if num is None:
        diff = position1 - position2
    else:
        diff = (position1 - position2)[:num]

    return diff.round(5)


if __name__ == "__main__":
    structure1_fn = sys.argv[1]
    structure2_fn = sys.argv[2]

    atoms1 = read(structure1_fn)
    atoms2 = read(structure2_fn)

    if len(sys.argv) == 3:
        diff = pos_diff(atoms1, atoms2)
    else:
        num = int(sys.argv[3])
        diff = pos_diff(atoms1, atoms2, num)

    print(
        f"Coordinate difference between {structure1_fn} with {structure2_fn}:"
    )
    print(diff)
