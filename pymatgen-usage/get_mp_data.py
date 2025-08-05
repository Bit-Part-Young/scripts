#!/usr/bin/env python3

"""获取 Materials Project 中的一元、二元数据"""

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


def get_mp_data(
    elements: str | list[str],
    stable: bool | None = None,
    energy_above_hull: float | None = None,
):
    """获取 Materials Project 中的一元、二元数据"""

    if isinstance(elements, list):
        # eg. ["Ti", "Al", "Ti-Al"]
        chemsys_list = elements + ["-".join(elements)]
    else:
        chemsys_list = elements

    if energy_above_hull is not None:
        energy_above_hull = (0.0, energy_above_hull)

    with MPRester(API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            chemsys=chemsys_list,
            fields=fields,
            is_stable=stable,
            energy_above_hull=energy_above_hull,
        )

    ndocs = len(docs)

    if isinstance(elements, list):
        print(f"\nTotally get {ndocs} {' '.join(elements)} MP binary data.\n")
    else:
        print(f"\nTotally get {ndocs} {elements} MP unary data.\n")

    data_list = []
    for doc in docs:

        data_dict = {
            "material_id": str(doc.material_id),
            "formula": doc.formula_pretty,
            "crystal_system": str(doc.symmetry.crystal_system),
            "spacegroup": doc.symmetry.symbol,
            "nsites": doc.nsites,
        }

        if isinstance(elements, list):
            structure_composition: dict[str, float] = doc.composition_reduced.as_dict()
            # 成分 分数形式
            structure_composition_frac = {
                element: structure_composition.get(element, 0.0)
                / sum(structure_composition.values())
                for element in elements
            }

            data_dict.update(structure_composition_frac)

        data_dict.update(
            {
                "fepa": doc.formation_energy_per_atom,
                "e_above_hull": doc.energy_above_hull,
                "stable": doc.is_stable,
            }
        )
        data_list.append(data_dict)

    df = pd.DataFrame(data_list).round(5)
    print(df)

    if isinstance(elements, list):
        csv_fn = f"{'_'.join(elements)}_mp.csv"
    else:
        csv_fn = f"{elements}_mp.csv"
    df.to_csv(csv_fn, index=False)

    print(f"\nData saved to {csv_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get Materials Project unary/binary data.",
        epilog="Author: SLY.",
    )

    parser.add_argument("elements", nargs="+", help="elements (e.g. Ti, Ti Al)")
    parser.add_argument("--stable", action="store_true", help="only get stable data")
    parser.add_argument(
        "-eah",
        "--energy_above_hull",
        type=float,
        metavar="energy_above_hull",
        help="maximum energy above hull",
    )

    args = parser.parse_args()

    if len(args.elements) == 1:
        args.elements = args.elements[0]
    elif len(args.elements) > 2:
        raise ValueError("Number of elements should not more than 2.")

    get_mp_data(
        elements=args.elements,
        stable=args.stable,
        energy_above_hull=args.energy_above_hull,
    )
