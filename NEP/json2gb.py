#!/usr/bin/env python3

"""将含 MP 晶界构型的 json 文件转换为 extxyz 格式"""

import argparse
import os

import pandas as pd
from ase.io import write
from monty.serialization import loadfn
from pymatgen.core.structure import Structure


def json2gb(json_fn: str):
    data = loadfn(json_fn)

    initial_output_fn = json_fn.replace(".json", "_initial.xyz")
    final_output_fn = json_fn.replace(".json", "_final.xyz")
    cif_output_fn = json_fn.replace(".json", "_cif.xyz")

    if os.path.exists(initial_output_fn):
        os.remove(initial_output_fn)
        os.remove(final_output_fn)
        os.remove(cif_output_fn)

    gb_info_list = []
    num_gb = 0
    for gb_data in data:
        # initial 和 final structure 有 grain_label property，具体为 top top_incident bottom bottom_incident
        # cif structure 无；cif 与 final structure 相同
        structure_initial: Structure = gb_data["initial_structure"]
        structure_final: Structure = gb_data["final_structure"]
        structure_cif = Structure.from_str(gb_data["cif"], fmt="cif")

        # 删除 grain_label 信息
        if "grain_label" in structure_initial.site_properties.keys():
            structure_initial.remove_site_property("grain_label")
            structure_final.remove_site_property("grain_label")

        natoms = len(structure_cif)
        gb_type = str(gb_data["type"])

        keys_list = [
            "sigma",
            "rotation_axis",
            "gb_plane",
            "rotation_angle",
            "gb_energy",
            "w_sep",
        ]
        gb_info = {key: gb_data[key] for key in keys_list}

        gb_info.update({"gb_type": gb_type})
        gb_info.update({"natoms": natoms})

        gb_plane = gb_data["gb_plane"]
        # 将 gb_plane 列表中的整数元素元素变为绝对值并从小到大拼接成字符串
        gb_plane_str = "".join([str(abs(index)) for index in sorted(gb_plane, key=abs)])
        gb_info.update({"gb_plane_str": gb_plane_str})

        gb_info_list.append(gb_info)

        write_atoms(structure_initial, gb_info, initial_output_fn)
        write_atoms(structure_final, gb_info, final_output_fn)
        write_atoms(structure_cif, gb_info, cif_output_fn)

        num_gb += 1

    gb_info_df = pd.DataFrame(gb_info_list).round(5)

    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    print(gb_info_df)

    csv_fn = json_fn.replace(".json", "_info.csv")
    gb_info_df.to_csv(csv_fn, index=False)

    print(
        f"Total {num_gb} GB configurations saved to {initial_output_fn} {final_output_fn} {cif_output_fn}."
    )


def write_atoms(structure: Structure, gb_info: dict, output_fn: str):
    atoms = structure.to_ase_atoms()
    atoms.info.update(gb_info)

    write(output_fn, atoms, format="extxyz", append=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert MP GB configurations from json to extxyz format.",
        epilog="Author: SLY.",
    )

    parser.add_argument("json_fn", help="json filename")
    args = parser.parse_args()

    json2gb(args.json_fn)
