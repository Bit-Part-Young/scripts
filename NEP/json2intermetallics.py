#!/usr/bin/env python3

"""将含 MP intermetallics 构型的 json 文件转换为 extxyz 格式，并保存相关构型数据"""

import argparse
import os

import pandas as pd
from ase.formula import Formula
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
        crystal_system = str(doc["symmetry"].crystal_system)
        is_stable = doc["is_stable"]
        material_id = doc["material_id"]
        formation_energy = doc["formation_energy_per_atom"]
        energy_above_hull = doc["energy_above_hull"]
        energy = doc["energy_per_atom"] * doc["nsites"]

        cell_info = atoms.cell.cellpar().round(3).tolist()
        natoms = len(atoms)
        formula = atoms.get_chemical_formula()

        composition = Formula(formula).count()
        element_list = list(composition.keys())
        # 成分，分数形式 {'Al': 0.4, 'Ti': 0.6}
        composition_fractional = {
            element: composition.get(element, 0) / natoms for element in element_list
        }

        info_dict = {
            "formula": formula,
            "natoms": natoms,
            "material_id": material_id,
            "energy": energy,
            "fe": formation_energy,
            "e_above_hull": energy_above_hull,
            "is_stable": is_stable,
            "crystal_system": crystal_system,
            "spg_number": spg_number,
            "spg_symbol": spg_symbol,
            "cell_info": cell_info,
        }

        info_dict.update(composition_fractional)

        info_list.append(info_dict)

        atoms_info_dict = info_dict.copy()
        # 删除不需要的 info
        keys_removed = ["cell_info", "energy"] + element_list
        for key in keys_removed:
            atoms_info_dict.pop(key, None)
        atoms_info_dict.update({"is_stable": str(is_stable)})
        atoms.info.update(atoms_info_dict)

        write(output_fn, atoms, format="extxyz", append=True)

    info_df = pd.DataFrame(info_list).round(5)
    print(info_df)

    csv_fn = json_fn.replace(".json", "_info.csv")
    info_df.to_csv(csv_fn, index=False)

    print("\nNote: The energy in MP is not real DFT calculated energy?")

    print(f"\n{output_fn} and {csv_fn} is generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert MP intermetallic configurations from json to extxyz format.",
        epilog="Author: SLY.",
    )

    parser.add_argument("json_fn", help="json filename")

    args = parser.parse_args()

    json2intermetallics(args.json_fn)
