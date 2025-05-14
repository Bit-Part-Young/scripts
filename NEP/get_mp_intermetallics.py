#!/usr/bin/env python3

"""从 MP 获取指定元素的 intermetallics 构型数据"""

import argparse
import os

from monty.serialization import dumpfn
from mp_api.client import MPRester

API_KEY = os.getenv("PMG_MAPI_KEY")


def get_mp_intermetallics(
    elements: list[str],
    formation_energy: bool = False,
):
    """
    从 MP 获取指定元素的 intermetallics 构型数据

    Args:
        elements: 元素列表
        formation_energy: 是否考虑形成能筛选条件
    """

    chemsys = "-".join(elements)

    fields = [
        "material_id",
        "is_stable",
        "formula_pretty",
        "formation_energy_per_atom",
        "energy_above_hull",
        "symmetry",
        "structure",
    ]

    kwargs = {}
    # 若 formation_energy 不为空，则search 函数添加 formation_energy 参数
    if formation_energy:
        formation_energy = (None, 0)
        kwargs["formation_energy"] = formation_energy

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

    print(f"Found {len(docs)} intermetallics for {chemsys}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get intermetallics from MP.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "-e",
        "--elements",
        nargs="+",
        type=str,
        help="elements",
    )

    parser.add_argument(
        "-fe",
        "--formation_energy",
        action="store_true",
        help="formation energy filter(less than 0.0 eV/atom)",
    )
    args = parser.parse_args()

    get_mp_intermetallics(args.elements, args.formation_energy)
