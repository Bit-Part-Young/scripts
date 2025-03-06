#!/usr/bin/env python3

"""获取 BCC/FCC/Diamond/HCP 晶体结构的最近邻距离"""

import argparse

import numpy as np
from ase.atoms import Atoms
from ase.build import bulk


def get_nn(crystalstructure: str) -> np.ndarray:
    """获取 BCC/FCC/Diamond/HCP 晶体结构的最近邻距离"""

    if crystalstructure in ["bcc", "fcc", "diamond"]:
        atoms = bulk("Mo", crystalstructure, a=1.0, cubic=True)
    elif crystalstructure == "hcp":
        atoms = bulk("Mo", crystalstructure, a=1.0)
    atoms_supercell: Atoms = atoms * (5, 5, 5)

    distances = atoms_supercell.get_all_distances(mic=True).round(5)
    nn_distances = np.unique(distances)
    nn_distances = nn_distances[distances > 0.0]

    print("  No.     R/a0        (R/a0)^2")
    for index, distance in enumerate(distances[:15], start=1):
        print(f"{index:>5}     {distance:<7}     {round(distance**2, 4):<7}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get the nearest neighbour distances for BCC/FCC/Diamond/HCP crystal structure.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "crystalstructure",
        type=str,
        choices=["bcc", "fcc", "diamond", "hcp"],
        help="crystal structure",
    )
    args = parser.parse_args()

    get_nn(args.crystalstructure)
