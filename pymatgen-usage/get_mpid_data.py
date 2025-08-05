#!/usr/bin/env python3

"""根据 Material Project material id 获取其数据"""

import argparse

import pandas as pd
from mp_api.client import MPRester
from pymatgen.core import SETTINGS


def get_mpid_data(material_id: str):
    """根据 Material Project material id 获取其数据"""

    PMG_API_KEY = SETTINGS.get("PMG_MAPI_KEY")

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

    with MPRester(PMG_API_KEY) as mpr:
        docs = mpr.materials.summary.search(fields=fields, material_ids=material_id)

    docs = docs[0]

    print(f"\n{material_id} data info:")
    data_dict = {
        "material_id": str(docs.material_id),
        "formula": docs.formula_pretty,
        "crystal_system": str(docs.symmetry.crystal_system),
        "spacegroup": docs.symmetry.symbol,
        "nsites": docs.nsites,
        "fepa": docs.formation_energy_per_atom,
        "energy_above_hull": docs.energy_above_hull,
        "is_stable": docs.is_stable,
    }

    df = pd.DataFrame([data_dict])
    print()
    print(df)

    return df


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get data by Material Project material id.",
        epilog="Author: SLY.",
    )

    parser.add_argument("mpid", help="Material Project material id (e.g. mp-149)")

    args = parser.parse_args()

    get_mpid_data(args.mpid)
