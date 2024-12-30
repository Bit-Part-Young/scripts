#!/usr/bin/env python3

"""获取构型中原子对的最小和最大距离"""

import argparse

import numpy as np
from ase.io import read


def get_distance(
    structure_fn: str = "POSCAR",
    mic: bool = False,
):
    """获取构型中原子对的最小和最大距离"""

    atoms = read(structure_fn)

    natoms = len(atoms)

    distances = atoms.get_all_distances(mic=mic)

    if natoms == 1:
        print("Number of atom is 1, no atomic pair.")

    elif natoms == 2 or len(np.unique(distances)) == 2:
        min_dist = np.unique(distances)[1]
        max_dist = min_dist

        print(f"Number of atom is {natoms} or all distances are same.")
        print(f"Distance: {round(min_dist, 5)} Å.")

    else:
        min_dist = np.unique(distances)[1]
        max_dist = np.unique(distances)[-1]

        print(f"Minimum distance: {round(min_dist, 5)} Å.")
        print(f"Maximum distance: {round(max_dist, 5)} Å.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get the minimum and maximum distance of atomic pair in a configuration.",
        epilog="Author: SLY",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        nargs="?",
        default="POSCAR",
        help="Structure filename",
    )

    parser.add_argument(
        "-mic",
        action="store_true",
        help="Use minimum image convention",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn
    mic = args.mic

    get_distance(
        structure_fn=args.structure_fn,
        mic=args.mic,
    )
