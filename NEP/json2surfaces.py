#!/usr/bin/env python3

"""从 json 文件中提取 MP 表面构型并保存为 xyz 文件"""

import argparse

import pandas as pd
from ase.io import write
from monty.serialization import loadfn
from pymatgen.core.structure import Structure


def json2surfaces(json_fn: str, max_index: int = None):
    data = loadfn(json_fn)

    surface_info_list = []
    num_surface = 0
    for surface in data[0]["surfaces"]:
        miller_index = surface["miller_index"]
        # 将负指数中的负号替换成 m 字符
        miller_index_str = "".join(
            map(lambda x: f"m{abs(x)}" if x < 0 else str(x), miller_index)
        )

        if max_index is None or abs(max(miller_index)) <= max_index:
            structure = Structure.from_str(surface["structure"], fmt="cif")
            atoms = structure.to_ase_atoms()

            keys_list = [
                "miller_index",
                "surface_energy",
                "surface_energy_EV_PER_ANG2",
                "is_reconstructed",
                "work_function",
                # "efermi",
                "area_fraction",
                # "has_wulff",
            ]
            surface_info = {key: surface[key] for key in keys_list}
            surface_info.update({"num_atoms": len(atoms)})
            surface_info.update({"miller_index": miller_index_str})
            surface_info_list.append(surface_info)

            atoms.info.update(surface_info)

            output_fn = json_fn.replace(".json", ".xyz")
            write(
                output_fn,
                atoms,
                format="extxyz",
                append=True,
            )

            num_surface += 1

    surface_info_df = pd.DataFrame(surface_info_list)

    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    print(surface_info_df)

    csv_fn = json_fn.replace(".json", "_surface_info.csv")
    surface_info_df.to_csv(csv_fn, index=False)

    print(f"\nTotal {num_surface} surface configurations saved to {output_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract MP surface configurations from json file to extxyz format.",
        epilog="Author: SLY.",
    )

    parser.add_argument("json_fn", type=str, help="json filename")
    parser.add_argument("-mi", "--max_index", type=int, help="max miller index")
    args = parser.parse_args()

    json2surfaces(args.json_fn, args.max_index)
