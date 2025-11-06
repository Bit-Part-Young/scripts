#!/usr/bin/env python3

"""从 MP 获取指定元素的 intermetallics 构型数据"""

import argparse
import os

from monty.serialization import dumpfn
from mp_api.client import MPRester

API_KEY = os.getenv("PMG_MAPI_KEY")


def get_mp_intermetallics(
    elements: list[str], formation_energy: bool = False, energy_above_hull: bool = False
):
    """
    从 MP 获取指定元素的 intermetallics 构型数据

    Args:
        elements: 元素列表
        formation_energy: 是否考虑形成能筛选条件
        energy_above_hull: 是否考虑能量高于 convex hull 筛选条件
    """

    chemsys = "-".join(elements)

    fields = [
        "material_id",
        "is_stable",
        "formula_pretty",
        "nsites",
        "energy_per_atom",
        "formation_energy_per_atom",
        "energy_above_hull",
        "symmetry",
        "structure",
    ]

    kwargs = {}
    # formation energy & energy above hull 的单位为 eV/atom
    if formation_energy:
        formation_energy = (None, 0)
        kwargs["formation_energy"] = formation_energy

    # energy above hull 上限取 0.04 eV/atom reference:
    # Efficient small-cell sampling for machine-learning potentials of multi-principal element alloys
    # http://arxiv.org/abs/2510.16697
    if energy_above_hull:
        energy_above_hull = (None, 0.04)
        kwargs["energy_above_hull"] = energy_above_hull

    with MPRester(api_key=API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            chemsys=chemsys,
            fields=fields,
            all_fields=False,
            **kwargs,
        )

    # 删除为请求数据为 null 的 key-value
    docs_simplified_list = []
    for doc in docs:
        doc_simplified = {}
        for field in fields:
            if getattr(doc, field) is not None:
                doc_simplified[field] = getattr(doc, field)

        docs_simplified_list.append(doc_simplified)

    json_fn = f"{chemsys}_mp_intermetallics.json"
    dumpfn(docs_simplified_list, json_fn)

    print("\nNote: The energy in MP is not real DFT calculated energy?")

    print(f"\nFound {len(docs)} intermetallics for {chemsys} and saved to {json_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get intermetallics from MP.", epilog="Author: SLY."
    )

    parser.add_argument("-e", "--elements", nargs="+", help="elements")

    parser.add_argument(
        "-fe",
        "--formation_energy",
        action="store_true",
        help="formation energy filter (less than 0.0 eV/atom)",
    )
    parser.add_argument(
        "-eah",
        "--energy_above_hull",
        action="store_true",
        help="energy above hull filter (less than 0.04 eV/atom)",
    )
    args = parser.parse_args()

    get_mp_intermetallics(args.elements, args.formation_energy, args.energy_above_hull)
