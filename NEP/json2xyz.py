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
    plusstress = json_data["outputs"]["virial_stress"]
    # 若 viral_stress 是 MTP cfg 提取的
    # virial = [
    #     plusstress[0],
    #     plusstress[5],
    #     plusstress[4],
    #     plusstress[5],
    #     plusstress[1],
    #     plusstress[3],
    #     plusstress[4],
    #     plusstress[3],
    #     plusstress[2],
    # ]
    # 若 viral_stress 是 VASP OUTCAR 提取的
    virial = [
        plusstress[0],
        plusstress[3],
        plusstress[5],
        plusstress[3],
        plusstress[1],
        plusstress[4],
        plusstress[5],
        plusstress[4],
        plusstress[2],
    ]
    atoms.info["virial"] = " ".join(map(str, virial))
    # 可能没有部分 key
    atoms.info["element"] = json_data.get("element", "False")
    atoms.info["num_atoms"] = json_data.get("num_atoms", "False")
    atoms.info["group"] = json_data.get("group", "False")
    atoms.info["description"] = json_data.get("description", "False")
    atoms.info["tag"] = json_data.get("tag", "False")

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
