#!/usr/bin/env python3

"""获取 BCC/FCC/Diamond 立方晶体结构的最近邻距离"""

import argparse

import numpy as np
import pandas as pd
from ase.build import bulk
from ase.geometry.geometry import get_distances


def get_nn(crystalstructure: str):
    """获取 BCC/FCC/Diamond 立方晶体结构的最近邻距离"""

    atoms = bulk("Mo", crystalstructure, a=1.0, cubic=True)
    atoms *= 10

    D, D_len = get_distances(
        np.diag(atoms.cell) / 2, atoms.positions, atoms.cell, atoms.pbc
    )

    center_index = D_len.argmin()
    center_position = atoms.positions[center_index]

    D, D_len = get_distances(center_position, atoms.positions, atoms.cell, atoms.pbc)
    D = D[0, :]
    D_len = D_len[0, :]

    decimals = 5
    D_len_rounded = np.round(D_len, decimals=decimals)
    unique_values, counts = np.unique(D_len_rounded, return_counts=True)

    df = pd.DataFrame(
        {
            "NNL": counts,
            "R/a0": unique_values,
            "(R/a0)^2": np.round(unique_values**2, 4),
        }
    ).round(5)
    df.index.name = "No."
    print(df.iloc[1:21, :])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get the nearest neighbour distances for BCC/FCC/Diamond cubic crystal structure.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "crystalstructure",
        choices=["bcc", "fcc", "diamond"],
        help="cubic crystal structure",
    )
    args = parser.parse_args()

    get_nn(args.crystalstructure)
