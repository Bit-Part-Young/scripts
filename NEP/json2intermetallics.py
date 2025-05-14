#!/usr/bin/env python3

"""从 json 文件中提取 MP intermetallics 构型并保存为 xyz 文件"""

import argparse
import os

import pandas as pd
from ase.io import write
from monty.serialization import loadfn
from pymatgen.core.structure import Structure


def json2intermetallics(json_fn: str):
    docs = loadfn(json_fn)

    output_fn = json_fn.replace(".json", ".xyz")
    if os.path.exists(output_fn):
        os.remove(output_fn)

    info_list = []
    for doc in docs:
        structure: Structure = doc["structure"]
        # 删除磁性信息
        if "magmom" in structure.site_properties.keys():
            structure.remove_site_property("magmom")
        atoms = structure.to_ase_atoms()

        spg_number = doc["symmetry"].number
        spg_symbol = doc["symmetry"].symbol
        is_stable = doc["is_stable"]
        material_id = doc["material_id"]
        formation_energy = doc["formation_energy_per_atom"]
        energy_above_hull = doc["energy_above_hull"]

        cell_info = atoms.cell.cellpar().round(3).tolist()
        natoms = len(atoms)
        formula = atoms.get_chemical_formula()

        info_dict = {
            "formula": formula,
            "natoms": natoms,
            "material_id": material_id,
            "fe": formation_energy,
            "e_above_hull": energy_above_hull,
            "is_stable": is_stable,
            "spg_number": spg_number,
            "spg_symbol": spg_symbol,
            "cell_info": cell_info,
        }

        info_list.append(info_dict)

        info_dict_new = info_dict.copy()
        info_dict_new.pop("cell_info")
        info_dict_new.update({"is_stable": str(is_stable)})
        atoms.info.update(info_dict_new)

        write(
            output_fn,
            atoms,
            format="extxyz",
            append=True,
        )

    info_df = pd.DataFrame(info_list).round(3)
    print(info_df)

    info_df.to_csv(json_fn.replace(".json", "_info.csv"), index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert MP intermetallics to xyz files.",
        epilog="Author: SLY.",
    )

    parser.add_argument("json_fn", type=str, help="json filename")

    args = parser.parse_args()

    json2intermetallics(args.json_fn)
