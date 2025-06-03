#!/usr/bin/env python3

"""使用 NEP 势函数进行 EOS 计算（晶格常数范围 ± 0.2 Å）"""

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


def eos_nep(
    symbol: str,
    a: float,
    covera: float = 1.633,
    potential_fn: str = "nep.txt",
):
    """使用 NEP 势函数 EOS 计算"""

    calc = CPUNEP(model_filename=potential_fn)

    scale = 0.2
    data = []
    for lc in np.arange(
        round(a - scale, 2),
        round(a + scale + 0.01, 2),
        0.01,
    ):
        lc = round(lc, 4)
        if symbol in ["Ti", "Zr"]:
            atoms = bulk(symbol, "hcp", a=lc, covera=covera)
        else:
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
    csv_fn = f"eos_nep_{symbol}.dat"
    df.to_csv(csv_fn, index=False, sep=" ")

    print(f"\nEOS data saved to {csv_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="EOS calculation using NEP potential.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Author: SLY.

Note: The lattice constant range ± 0.2 Å, interval 0.01 Å, total 41 data points.

Example:
    python eos_nep.py Ti 2.9 1.582
    python eos_nep.py Nb 3.3
        """,
    )

    parser.add_argument("symbol", type=str, help="element symbol")
    parser.add_argument("a", type=float, help="a (rough value)")

    parser.add_argument(
        "covera",
        type=float,
        default=1.633,
        const=1.633,
        nargs="?",
        help="c/a (default: 1.633)",
    )

    args = parser.parse_args()

    eos_nep(args.symbol, args.a, args.covera)
