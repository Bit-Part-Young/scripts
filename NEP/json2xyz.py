"""
将 json 构型及其数据文件转换为 xyz 文件

文件格式: [{}, {}, ...]

dict_keys(['element', 'num_atoms', 'group', 'description', 'tag', 'structure', 'outputs'])

outputs 数据中的 key
dict_keys(['energy', 'forces', 'virial_stress'])
"""

import numpy as np
from ase.io import write
from monty.serialization import loadfn
from pymatgen.core.structure import Structure


def json2xyz(json_data: dict):
    """将 json 结构文件转换为 xyz 文件"""

    structure: Structure = json_data["structure"]
    atoms = structure.to_ase_atoms()

    atoms.arrays["forces"] = np.array(json_data["outputs"]["forces"])
    atoms.info["energy"] = json_data["outputs"]["energy"]

    # 将 virial_stress 转换为 virial；需确定分量顺序
    virial_stress = json_data["outputs"]["virial_stress"]
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
    keys = ["element", "num_atoms", "group", "description", "tag"]
    for key in keys:
        atoms.info[key] = json_data.get(key, "False")

    return atoms


def main():

    json_fn = "example.json"
    xyz_fn = "example.xyz"

    json_data_list = loadfn(json_fn)

    flag = 0
    for json_data in json_data_list:
        atoms = json2xyz(json_data)

        write(xyz_fn, atoms, format="extxyz", append=True)

        flag += 1

        print(f"No. {flag} structure has been processed.")


if __name__ == "__main__":
    main()
