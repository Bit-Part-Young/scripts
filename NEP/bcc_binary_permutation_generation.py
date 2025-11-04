"""
生成 BCC 2x2x2 超胞的二元固溶体的不同浓度（1~15/16）的对称性非等同的置换构型
每种二元固溶体共 329 个置换构型
核心工具：bsym Python 包
reference: https://www.nature.com/articles/s41524-024-01330-6

"""

import os
from itertools import combinations

from ase.build import bulk
from ase.io import write
from bsym.interface.pymatgen import unique_structure_substitutions
from pymatgen.core.structure import Structure

element_lc_dict = {
    "Ti": 3.252,
    "Al": 3.232,
    "Nb": 3.307,
    "Mo": 3.164,
    "Zr": 3.572,
    "V": 2.996,
}


def get_binary_lc(natoms_list: list, element_list: list):
    """获取 BCC 2x2x2 超胞的二元合金不同浓度（1~15/16）的晶格常数"""

    lc = (
        element_lc_dict[element_list[0]] * natoms_list[0]
        + element_lc_dict[element_list[1]] * natoms_list[1]
    ) / sum(natoms_list)

    return round(lc, 3)


def permutation_generation(natoms_list: list, element_list: list):
    """生成置换对称性非等同的构型"""

    lc = get_binary_lc(natoms_list, element_list)
    atoms_initial = bulk("Mo", "bcc", a=lc, cubic=True) * (2, 2, 2)
    # atoms_initial = bulk("Mo", "bcc", a=3.254, cubic=True) * (2, 2, 2)
    structure_initial = Structure.from_ase_atoms(atoms_initial)

    subs_structures = unique_structure_substitutions(
        structure=structure_initial,
        to_substitute="Mo",
        site_distribution={
            element_list[0]: natoms_list[0],
            element_list[1]: natoms_list[1],
        },
        verbose=True,
        show_progress=True,
    )

    return subs_structures


if __name__ == "__main__":

    elements = ["Ti", "Al", "Nb", "Mo", "Zr", "V"]

    xyz_total_fn = "bcc_binary_permutation_binary.xyz"
    if os.path.exists(xyz_total_fn):
        os.remove(xyz_total_fn)

    flag = 0
    for element1, element2 in combinations(elements, 2):

        chemsys = "-".join([str(flag)] + [element1, element2])
        element_list = [element1, element2]

        os.makedirs(chemsys, exist_ok=True)

        os.chdir(chemsys)

        print(f"Binary {chemsys}:")
        for i in range(1, 16):
            os.makedirs(f"{i}", exist_ok=True)
            os.chdir(f"{i}")

            natoms_list = [i, 16 - i]
            subs_structures = permutation_generation(natoms_list, element_list)
            nstructures = len(subs_structures)

            xyz_fn = f"{element1}{element2}_binary_permutation_conc_{i}.xyz"
            if os.path.exists(xyz_fn):
                os.remove(xyz_fn)

            for subs_structure in subs_structures:
                cwd_path = os.getcwd()
                atoms = Structure.to_ase_atoms(subs_structure)

                write(xyz_fn, atoms, format="extxyz", append=True)

                xyz_total_fn = f"../../{xyz_total_fn}"
                write(xyz_total_fn, atoms, format="extxyz", append=True)

            os.chdir("..")

            print(f"concentration {i}/16 done.")

        os.chdir("..")
