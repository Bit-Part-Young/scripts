#!/usr/bin/env python3

"""
生成 FCC 1x2x2 超胞（16 个原子）的二元固溶体的不同浓度（0.125~0.875 共 7 个浓度）的对称性非等同的置换构型
每种 BCC 二元固溶体中每种浓度至多保留 20 个构型，共 110 个置换构型

FCC 超胞选择 1x2x2 同 reference 中的做法

需要编译 enumlib，并将 enum.x 和 aux_src/makeStr.py 添加到 PATH 中

reference: http://link.aps.org/supplemental/10.1103/PhysRevMaterials.7.103602
"""

import os
import subprocess
from itertools import combinations

import numpy as np
from ase.build import bulk
from ase.io import write
from pymatgen.core.structure import Structure
from pymatgen.transformations.advanced_transformations import (
    EnumerateStructureTransformation,
)


def get_fcc_binary_lc(natoms_list: list, element_list: list):
    # [ ] FCC 1x2x2 这样计算是否合适？
    """获取 FCC 1x2x2 超胞的二元合金不同浓度的平均晶格常数"""

    lc = (
        element_lc_dict[element_list[0]] * natoms_list[0]
        + element_lc_dict[element_list[1]] * natoms_list[1]
    ) / sum(natoms_list)

    return round(lc, 3)


def enumerate_solution(element_initial: str, element_replace: str, num_kept: int = 20):

    est = EnumerateStructureTransformation()

    natoms = len(bulk("Al", "fcc", a=4.041, cubic=True) * (1, 2, 2))
    # FCC 0.125~0.875 共 7 个浓度
    concentrations = np.arange(0.125, 1, 0.125)
    num_replace_list = [int(natoms * conc) for conc in concentrations]

    alloy_structures = []
    for i, num_replace in enumerate(num_replace_list):
        natoms_list = [natoms - num_replace, num_replace]
        lc = get_fcc_binary_lc(natoms_list, [element_initial, element_replace])

        atoms = bulk(element_initial, "fcc", a=lc, cubic=True)
        atoms = atoms * (1, 2, 2)
        initial_structure = Structure.from_ase_atoms(atoms)

        os.makedirs(f"{i+1}", exist_ok=True)
        os.chdir(f"{i+1}")

        xyz_fn = (
            f"{element_initial}{element_replace}_fcc_binary_permutation_conc_{i+1}.xyz"
        )
        if os.path.exists(xyz_fn):
            os.remove(xyz_fn)

        s_copy = initial_structure.copy()
        s_copy.replace_species(
            {
                element_initial: {
                    element_replace: num_replace / natoms,
                    element_initial: (natoms - num_replace) / natoms,
                }
            }
        )

        # 返回前 100 个低能量构型；不添加 return_ranked_list 则返回能量最低的构型；类型为 list[dict]
        ss = est.apply_transformation(s_copy, return_ranked_list=100)
        ss = [s["structure"] for s in ss]
        # 至多保留 N 个构型
        ss_kept = ss[:num_kept]

        alloy_structures.append(ss_kept)

        atoms_list = [s.to_ase_atoms() for s in ss_kept]
        write(xyz_fn, atoms_list)

        print(
            "take {} from total {} of {}.".format(
                len(alloy_structures[i]), len(ss), ss[0].formula
            )
        )

        print(f"conc {concentrations[i]} done.\n")

        os.chdir("..")

    return alloy_structures


if __name__ == "__main__":

    # FCC 结构元素的晶格常数
    element_lc_dict = {
        "Ti": 4.108,
        "Al": 4.041,
        "Nb": 4.215,
        "Mo": 4.004,
        "Zr": 4.530,
        "V": 3.817,
    }

    elements = ["Ti", "Al", "Nb", "Mo", "Zr", "V"]

    flag = 0
    for element1, element2 in combinations(elements, 2):
        flag += 1

        chemsys = "-".join([str(flag)] + [element1, element2])

        os.makedirs(chemsys, exist_ok=True)

        os.chdir(chemsys)

        print(f"FCC Binary {chemsys}:")

        xyz_fn = f"{element1}{element2}_fcc_binary_permutation.xyz"
        if os.path.exists(xyz_fn):
            os.remove(xyz_fn)

        enumerate_solution(element1, element2, num_kept=20)

        subprocess.run(f"cat */*.xyz > {xyz_fn}", shell=True)

        os.chdir("..")
