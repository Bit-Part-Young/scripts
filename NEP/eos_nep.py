#!/usr/bin/env python3

"""使用 NEP 势函数 EOS 计算"""

import argparse

import numpy as np
import pandas as pd
from ase.build import bulk
from calorine.calculators import CPUNEP, GPUNEP

crystalstructure_dict = {
    "Al": "fcc",
    "Nb": "bcc",
    "Mo": "bcc",
    "V": "bcc",
    "Zr": "hcp",
    "Ti": "hcp",
}


# [ ]  添加 HCP 结构 EOS
def eos_nep(symbol: str, lc_round: float, potential_fn: str = "nep.txt"):
    """使用 NEP 势函数 EOS 计算"""

    calc = CPUNEP(model_filename=potential_fn)

    scale = 0.2
    data = []
    for lc in np.arange(
        round(lc_round - scale, 2),
        round(lc_round + scale, 2),
        0.01,
    ):
        lc = round(lc, 4)
        atoms = bulk(symbol, crystalstructure_dict[symbol], a=lc, cubic=True)

        atoms.calc = calc

        natoms = len(atoms)
        energy_pa = round(atoms.get_potential_energy() / natoms, 5)
        volume_pa = round(atoms.get_volume() / natoms, 5)

        data_dict = {
            "lc": lc,
            "volume_pa": volume_pa,
            "energy_pa": energy_pa,
        }
        data.append(data_dict)

    df = pd.DataFrame(data)
    df.to_csv(f"ev_{symbol}_nep.dat", index=False, sep=" ")

    print(f"EOS data saved to ev_{symbol}.dat.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="EOS calculation using NEP potential.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("symbol", type=str, help="element symbol")
    parser.add_argument("lc_round", type=float, help="rough lattice constant")

    args = parser.parse_args()

    eos_nep(args.symbol, args.lc_round)
