#!/usr/bin/env python3

"""获取构型中原子对的最小和最大距离"""

import argparse

import numpy as np
from ase.atoms import Atoms
from ase.io import read


def get_distance(atoms: Atoms, mic=False):
    """获取构型中原子对的最小和最大距离"""

    distances = atoms.get_all_distances(mic=mic).round(5)

    min_dist = np.min(distances[distances > 0.0])
    max_dist = np.max(distances)

    return min_dist, max_dist


def main(
    structure_fn: str = "POSCAR",
    mic: bool = False,
):

    fn_suffix = structure_fn.split(".")[-1]
    if fn_suffix == "xyz":
        atoms_list = read(structure_fn, index=":", format="extxyz")

        min_dist_list, max_dist_list = [], []
        for atoms in atoms_list:
            if len(atoms) == 1:
                print("Number of atom is 1, pass.")
            else:
                min_dist, max_dist = get_distance(atoms, mic=mic)

            min_dist_list.append(min_dist)
            max_dist_list.append(max_dist)

        min_dist_global = np.min(min_dist_list)
        max_dist_global = np.max(max_dist_list)

        print(f"Global Minimum distance: {min_dist_global} Å.")
        print(f"Global Maximum distance: {max_dist_global} Å.")

    else:
        atoms = read(structure_fn)

        natoms = len(atoms)

        if natoms == 1:
            print("Number of atom is 1, no atomic pair, pass.")
        else:
            min_dist, max_dist = get_distance(atoms, mic=mic)

            print(f"Minimum distance: {min_dist} Å.")
            print(f"Maximum distance: {max_dist} Å.")

            if natoms == 2 or np.isclose(min_dist, max_dist):
                print(f"Number of atom is {natoms} or all distances are same.")


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

    main(
        structure_fn=args.structure_fn,
        mic=args.mic,
    )
