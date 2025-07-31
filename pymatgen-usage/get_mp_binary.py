#!/usr/bin/env python3

"""获取 Materials Project 二元数据"""

import argparse
import os

import pandas as pd
from mp_api.client import MPRester

fields = [
    "material_id",
    "structure",
    "composition_reduced",
    "formula_pretty",
    "symmetry",
    "nsites",
    "energy_per_atom",
    "formation_energy_per_atom",
    "energy_above_hull",
    "is_stable",
]

API_KEY = os.getenv("PMG_MAPI_KEY")


def get_mp_binary(
    element_list: list[str],
    stable: bool = False,
    energy_above_hull: float | None = None,
):
    """获取 Materials Project 二元数据"""

    if len(element_list) != 2:
        raise ValueError("element_list must be a list of two elements.")

    chemsys_list = element_list + ["-".join(element_list)]

    with MPRester(API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            chemsys=chemsys_list,  # eg. ["Ti", "Al", "Ti-Al"]
            fields=fields,
            is_stable=stable,
            energy_above_hull=(0.0, energy_above_hull),
        )

    ndocs = len(docs)
    print(f"\nTotally get {ndocs} {' '.join(element_list)} MP binary data.\n")

    data_list = []
    for doc in docs:
        structure_composition = doc.composition_reduced.as_dict()
        # 成分 分数形式
        structure_composition_frac = {
            element: structure_composition.get(element, 0)
            / sum(structure_composition.values())
            for element in element_list
        }

        data_dict = {
            "material_id": str(doc.material_id),
            "formula": doc.formula_pretty,
            "crystal_system": str(doc.symmetry.crystal_system),
            "spacegroup": doc.symmetry.symbol,
            "nsites": doc.nsites,
        }
        data_dict.update(structure_composition_frac)
        data_dict.update(
            {
                # "epa": doc.energy_per_atom,  # 非平均原子能量
                "fepa": doc.formation_energy_per_atom,
                "e_above_hull": doc.energy_above_hull,
                "stable": doc.is_stable,
            }
        )
        data_list.append(data_dict)

    df = pd.DataFrame(data_list).round(5)
    print(df)

    csv_fn = f"{'_'.join(element_list)}.csv"
    df.to_csv(csv_fn, index=False)

    print(f"\nData saved to {csv_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get Materials Project binary data.",
        epilog="Author: SLY.",
    )

    parser.add_argument("element_list", nargs=2, help="element list (e.g. Ti Al)")
    parser.add_argument("--stable", action="store_true", help="only get stable data")
    parser.add_argument(
        "-eah",
        "--energy_above_hull",
        type=float,
        metavar="energy_above_hull",
        help="maximum energy above hull",
    )

    args = parser.parse_args()

    get_mp_binary(
        element_list=args.element_list,
        stable=args.stable,
        energy_above_hull=args.energy_above_hull,
    )
