#!/usr/bin/env python3

"""
将 json 构型及其数据文件转换为 xyz 文件

文件格式: [{}, {}, ...]

dict_keys(['element', 'num_atoms', 'group', 'description', 'tag', 'structure', 'outputs'])

outputs 数据中的 key
dict_keys(['energy', 'forces', 'virial_stress'])
"""

import argparse
import os
from typing import Any

import numpy as np
from ase.io import write
from monty.serialization import loadfn
from pymatgen.core.structure import Structure


def json2atoms(json_data: dict, output_virial: bool = False):
    """将 json 结构文件转换为 xyz 文件"""

    structure: Structure = json_data["structure"]
    atoms = structure.to_ase_atoms()

    outputs: dict[str, Any] = json_data.get("outputs", {})
    if outputs is not None:
        if "forces" in outputs.keys():
            atoms.arrays["forces"] = np.array(outputs.get("forces", []))
        if "energy" in outputs.keys():
            atoms.info["energy"] = outputs.get("energy", 0.0)
        if "virial_stress" in outputs.keys() and output_virial:
            # 将 virial_stress 转换为 virial；需确定分量顺序
            virial_stress = outputs.get("virial_stress", [])
            # Voigt order
            virial = [
                virial_stress[0],
                virial_stress[5],
                virial_stress[4],
                virial_stress[5],
                virial_stress[1],
                virial_stress[3],
                virial_stress[4],
                virial_stress[3],
                virial_stress[2],
            ]
            # 非 Voigt order
            virial = [
                virial_stress[0],
                virial_stress[3],
                virial_stress[5],
                virial_stress[3],
                virial_stress[1],
                virial_stress[4],
                virial_stress[5],
                virial_stress[4],
                virial_stress[2],
            ]
            atoms.info["virial"] = " ".join(map(str, virial))

    # 可能没有部分 key
    additional_keys = ["element", "num_atoms", "group", "description", "tag"]
    for key in additional_keys:
        atoms.info[key] = json_data.get(key, "False")

    return atoms


def json2xyz(json_fn: str, xyz_fn: str, output_virial: bool = False):
    json_data_list = loadfn(json_fn)

    if os.path.exists(xyz_fn):
        os.remove(xyz_fn)

    flag = 0
    for json_data in json_data_list:
        atoms = json2atoms(json_data, output_virial)

        write(xyz_fn, atoms, format="extxyz", append=True)

        flag += 1

        print(f"No. {flag} structure has been processed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert json to extxyz format.")

    parser.add_argument("json_fn", help="inputjson filename")
    parser.add_argument("xyz_fn", help="output extxyz filename")
    parser.add_argument(
        "-virial", action="store_true", help="whether to output virial / stress"
    )
    args = parser.parse_args()

    json2xyz(args.json_fn, args.xyz_fn, args.virial)
