#!/usr/bin/env python3

"""获取 表面性质 和 晶界 数据"""

import argparse
import os

from monty.serialization import dumpfn
from mp_api.client import MPRester

API_KEY = os.getenv("PMG_MAPI_KEY")


def get_surface_properties(mp_id: str, system: str):
    """获取 表面性质 数据"""

    with MPRester(api_key=API_KEY) as mpr:
        surface_properties_docs = mpr.materials.surface_properties.search(
            material_ids=mp_id
        )

    num_surface = len(surface_properties_docs[0].surfaces)
    print(f"{mp_id} {system} has {num_surface} surfaces data.\n")

    json_fn = f"{mp_id}_{system}_surface.json"
    dumpfn(surface_properties_docs, json_fn)

    print(f"Surface properties data saved to {json_fn}.")


def get_grain_boundary(mp_id: str, system: str):
    """获取 晶界 数据"""

    with MPRester(api_key=API_KEY) as mpr:
        grain_boundary_docs = mpr.materials.grain_boundaries.search(material_ids=mp_id)

    num_gb = len(grain_boundary_docs)
    print(f"{mp_id} {system} has {num_gb} grain boundary data.\n")

    json_fn = f"{mp_id}_{system}_gb.json"
    dumpfn(grain_boundary_docs, json_fn)

    print(f"Grain boundary data saved to {json_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get surface properties and grain boundary data from MP with given material id and system.",
        epilog="Author: SLY.",
    )

    parser.add_argument("mp_id", type=str, help="MP id")
    parser.add_argument("system", type=str, help="System element symbol")
    args = parser.parse_args()

    get_surface_properties(args.mp_id, args.system)
    print("\n" * 2)
    get_grain_boundary(args.mp_id, args.system)
