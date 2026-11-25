#!/usr/bin/env python3

"""
获取 BCC/FCC/Diamond/HCP 晶体结构的各近邻壳层距离（以晶格常数 a 为单位）及配位数

FCC 近邻: 3rd 24, 7th 48, 13th 72
HCP: 默认理想轴比 c/a = sqrt(8/3)；非理想轴比时第 1 近邻壳层劈裂为 6 + 6

reference: https://www.ifmpan.poznan.pl/~urbaniak/fcc%20and%20its%20neigbours01.pdf
"""

import argparse

import numpy as np
import pandas as pd
from ase.build import bulk
from ase.neighborlist import neighbor_list


def get_nn(
    crystalstructure: str, nshells: int = 20, covera: float = np.sqrt(8 / 3), a: float | None = None
) -> pd.DataFrame:
    """获取晶体结构前 nshells 个近邻壳层的距离 R/a 及各壳层原子数；给定 a 时额外输出 R (Å)"""

    kwargs = {"covera": covera} if crystalstructure == "hcp" else {"cubic": True}
    atoms = bulk("X", crystalstructure, a=1.0, **kwargs)

    # 逐步增大截断半径直至覆盖 nshells 个壳层（各原子等价，取 0 号原子为中心）
    cutoff = 2.0
    while True:
        i, d = neighbor_list("id", atoms, cutoff)
        shells, counts = np.unique(np.round(d[i == 0], 5), return_counts=True)
        if len(shells) > nshells:  # 最外层可能不完整，多取一层再截断
            break
        cutoff += 1.0

    df = pd.DataFrame(
        {"count": counts, "R/a": shells, "(R/a)^2": np.round(shells**2, 4)},
        index=pd.RangeIndex(1, len(shells) + 1, name="shell"),
    ).iloc[:nshells]
    if a is not None:
        df["R(Å)"] = np.round(df["R/a"] * a, 4)

    print(df.to_string())

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get the neighbour shell distances and counts of BCC/FCC/Diamond/HCP crystals.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument(
        "crystalstructure", choices=["bcc", "fcc", "diamond", "hcp"], help="crystal structure"
    )
    parser.add_argument("-n", "--nshells", type=int, default=20, help="number of neighbour shells")
    parser.add_argument(
        "-a", type=float, metavar="A", help="lattice constant a (Å) to also print R in Å"
    )
    parser.add_argument(
        "--covera", type=float, default=np.sqrt(8 / 3), help="c/a ratio for hcp (ideal: sqrt(8/3))"
    )
    args = parser.parse_args()

    get_nn(args.crystalstructure, nshells=args.nshells, covera=args.covera, a=args.a)
