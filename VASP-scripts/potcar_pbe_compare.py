"""比较 VASP 和 pymatgen 推荐的常用元素 PBE 赝势 (VASP5.4)"""

import sys

# Data source: https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/io/vasp/MPRelaxSet.yaml
pymatgen_pbe_dict = {
    "H": "H",
    "He": "He",
    "Li": "Li_sv",
    "Be": "Be_sv",
    "B": "B",
    "C": "C",
    "N": "N",
    "O": "O",
    "Mg": "Mg_pv",
    "Al": "Al",
    "Si": "Si",
    "S": "S",
    "Cl": "Cl",
    "Ca": "Ca_sv",
    "Sc": "Sc_sv",
    "Ti": "Ti_pv",
    "V": "V_pv",
    "Cr": "Cr_pv",
    "Mn": "Mn_pv",
    "Fe": "Fe_pv",
    "Co": "Co",
    "Ni": "Ni_pv",
    "Cu": "Cu_pv",
    "Zn": "Zn",
    "Zr": "Zr_sv",
    "Nb": "Nb_pv",
    "Mo": "Mo_pv",
    "Sn": "Sn_d",
    "Hf": "Hf_pv",
    "Ta": "Ta_pv",
    "W": "W_sv",
}

# Data source: https://www.vasp.at/wiki/index.php/Available_PAW_potentials
vasp_pbe_dict = {
    "H": "H",
    "He": "He",
    "Li": "Li_sv",
    "Be": "Be",
    "B": "B",
    "C": "C",
    "N": "N",
    "O": "O",
    "Mg": "Mg",
    "Al": "Al",
    "Si": "Si",
    "S": "S",
    "Cl": "Cl",
    "Ca": "Ca_sv",
    "Sc": "Sc_sv",
    "Ti": "Ti_sv",
    "V": "V",
    "Cr": "Cr_pv",
    "Mn": "Mn_pv",
    "Fe": "Fe",
    "Co": "Co",
    "Ni": "Ni",
    "Cu": "Cu",
    "Zn": "Zn",
    "Zr": "Zr_sv",
    "Nb": "Nb_sv",
    "Mo": "Mo_sv",
    "Sn": "Sn_d",
    "Hf": "Hf_pv",
    "Ta": "Ta_pv",
    "W": "W_sv",
}


def potcar_pbe_compare(element: str):
    """比较 VASP 和 pymatgen 推荐的常用元素 PBE 赝势(VASP5.4)"""

    print(f"VASP recommended: {vasp_pbe_dict[element]}")
    print(f"Pymatgen recommended: {pymatgen_pbe_dict[element]}")


if __name__ == "__main__":

    element = sys.argv[1]
    potcar_pbe_compare(element=element)
