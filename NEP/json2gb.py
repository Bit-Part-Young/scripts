#!/usr/bin/env python3

"""从 json 文件中提取 MP 晶界构型并保存为 xyz 文件"""

import argparse

import pandas as pd
from ase.io import write
from monty.serialization import loadfn
from pymatgen.core.structure import Structure


def json2gb(json_fn: str):
    data = loadfn(json_fn)

    gb_info_list = []
    num_gb = 0
    for gb_data in data:
        initial_structure = gb_data["initial_structure"]
        final_structure = gb_data["final_structure"]
        structure_cif = Structure.from_str(gb_data["cif"], fmt="cif")

        keys_list = [
            "sigma",
            "rotation_axis",
            "gb_plane",
            "rotation_angle",
            "gb_energy",
            "w_sep",
        ]
        gb_info = {key: gb_data[key] for key in keys_list}

        gb_info.update({"gb_type": str(gb_data["type"])})
        gb_info.update({"natoms": len(structure_cif)})

        gb_info_list.append(gb_info)

        initial_output_fn = json_fn.replace(".json", "_initial.xyz")
        final_output_fn = json_fn.replace(".json", "_final.xyz")
        cif_output_fn = json_fn.replace(".json", "_cif.xyz")

        write_atoms(initial_structure, gb_info, initial_output_fn)
        write_atoms(final_structure, gb_info, final_output_fn)
        write_atoms(structure_cif, gb_info, cif_output_fn)

        num_gb += 1

    gb_info_df = pd.DataFrame(gb_info_list)

    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    print(gb_info_df)

    csv_fn = json_fn.replace(".json", "_gb_info.csv")
    gb_info_df.to_csv(csv_fn, index=False)

    print(
        f"Total {num_gb} surface configurations saved to {initial_output_fn} {final_output_fn} {cif_output_fn}."
    )


def write_atoms(structure: Structure, gb_info: dict, output_fn: str):
    atoms = structure.to_ase_atoms()
    atoms.info.update(gb_info)

    write(
        output_fn,
        atoms,
        format="extxyz",
        append=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract MP surface configurations from json file to extxyz format.",
        epilog="Author: SLY.",
    )

    parser.add_argument("json_fn", type=str, help="json filename")
    args = parser.parse_args()

    json2gb(args.json_fn)
