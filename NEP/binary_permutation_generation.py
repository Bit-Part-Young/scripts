#!/usr/bin/env python3

"""
生成 BCC/FCC/HCP 2x2x2 超胞的二元固溶体的不同浓度的对称性非等同的置换构型
需要编译 enumlib，并将 enum.x 和 aux_src/makeStr.py 添加到 PATH 中

reference: http://link.aps.org/supplemental/10.1103/PhysRevMaterials.7.103602
"""

import argparse

import numpy as np
from ase.build import bulk
from pymatgen.core.structure import Structure
from pymatgen.transformations.advanced_transformations import (
    EnumerateStructureTransformation,
)


def enumerate_solution(initial_structure: Structure):
    est = EnumerateStructureTransformation()

    alloy_structures = []

    natoms = initial_structure.num_sites

    print(
        initial_structure.formula,
        initial_structure.lattice.abc,
        initial_structure.lattice.angles,
    )

    # BCC
    concentrations = np.arange(0.0625, 1, 0.0625)
    # FCC
    # ratios = np.arange(0.125, 1, 0.125)
    num_replace_list = [int(natoms * conc) for conc in concentrations]

    # [ ] 替换元素部分待优化
    element_initial = str(initial_structure.composition.elements[0])
    if element_initial == "Al":
        element_replace = "Ti"
    elif element_initial == "Ti":
        element_replace = "Al"

    for i, num_replace in enumerate(num_replace_list):
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

        # [ ] 这里需要优化，只保留前 30 个构型
        alloy_structures.append(ss[:30])

        print(
            "take {} from total {} of {}".format(
                len(alloy_structures[i]), len(ss), ss[0].formula
            )
        )

    return alloy_structures


if __name__ == "__main__":
    # atoms = bulk("Al", "fcc", a=4.041, cubic=True)
    atoms = bulk("Al", "bcc", a=3.307, cubic=True)
    atoms = atoms * (2, 2, 2)

    initial_structure = Structure.from_ase_atoms(atoms)
    enumerate_solution(initial_structure)
