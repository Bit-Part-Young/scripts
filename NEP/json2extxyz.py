"""
将 json 构型及其数据文件转换为 extxyz 格式

文件格式: [{}, {}, ...]

dict_keys(['element', 'num_atoms', 'group', 'description', 'tag', 'structure', 'outputs'])

outputs 数据中的 key
dict_keys(['energy', 'forces', 'virial_stress'])
"""

import json

import numpy as np
from ase.io import write
from pymatgen.core.structure import Structure

# json_fn = "example.json"
# extxyz_fn = "example.xyz"
json_fn = "quinary.json"
extxyz_fn = "quinary.xyz"


with open(json_fn, "r") as f:
    json_data_list = json.load(f)


def json2extxyz(json_data: dict):
    """将 json 结构文件转换为 extxyz 格式"""

    structure = Structure.from_dict(json_data["structure"])
    atoms = structure.to_ase_atoms()

    atoms.arrays["forces"] = np.array(json_data["outputs"]["forces"])
    atoms.info["energy"] = json_data["outputs"]["energy"]
    atoms.info["virial_stress"] = " ".join(map(str, json_data["outputs"]["virial_stress"]))
    # 可能没有部分 key
    atoms.info["element"] = json_data.get("element", "False")
    atoms.info["num_atoms"] = json_data.get("num_atoms", "False")
    atoms.info["group"] = json_data.get("group", "False")
    atoms.info["description"] = json_data.get("description", "False")
    atoms.info["tag"] = json_data.get("tag", "False")

    return atoms


def main():
    flag = 0
    for json_data in json_data_list:
        atoms = json2extxyz(json_data)

        write(
            extxyz_fn,
            atoms,
            format="extxyz",
            append=True,
        )

        flag += 1

        print(f"No. {flag} structure has been processed.")


if __name__ == "__main__":
    main()
