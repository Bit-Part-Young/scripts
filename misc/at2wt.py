"""将 Ti-22Al-23Nb-1Mo-1Zr 格式化学式原子百分比转化成质量百分比"""

import re
from typing import Dict

# atomic mass data: got from pymatgen
atomic_mass_dict = {
    "H": 1.00794,
    "He": 4.002602,
    "Li": 6.941,
    "Be": 9.012182,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984032,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.305,
    "Al": 26.9815386,
    "Si": 28.0855,
    "P": 30.973762,
    "S": 32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955912,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938045,
    "Fe": 55.845,
    "Co": 58.933195,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.409,
    "Ga": 69.723,
    "Ge": 72.64,
    "As": 74.9216,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90585,
    "Zr": 91.224,
    "Nb": 92.90638,
    "Mo": 95.94,
    "Tc": 98.0,
    "Ru": 101.07,
    "Rh": 102.9055,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.71,
    "Sb": 121.76,
    "Te": 127.6,
    "I": 126.90447,
    "Xe": 131.293,
    "Cs": 132.9054519,
    "Ba": 137.327,
    "La": 138.90547,
    "Ce": 140.116,
    "Pr": 140.90765,
    "Nd": 144.242,
    "Pm": 145.0,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.92535,
    "Dy": 162.5,
    "Ho": 164.93032,
    "Er": 167.259,
    "Tm": 168.93421,
    "Yb": 173.04,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.94788,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.966569,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.9804,
    "Po": 210.0,
    "At": 210.0,
    "Rn": 220.0,
    "Fr": 223.0,
    "Ra": 226.0,
    "Ac": 227.0,
    "Th": 232.03806,
    "Pa": 231.03588,
    "U": 238.02891,
    "Np": 237.0,
    "Pu": 244.0,
    "Am": 243.0,
    "Cm": 247.0,
    "Bk": 247.0,
    "Cf": 251.0,
    "Es": 252.0,
    "Fm": 257.0,
    "Md": 258.0,
    "No": 259.0,
    "Lr": 262.0,
    "Rf": 267.0,
    "Db": 268.0,
    "Sg": 269.0,
    "Bh": 270.0,
    "Hs": 270.0,
    "Mt": 278.0,
    "Ds": 281.0,
    "Rg": 282.0,
    "Cn": 285.0,
    "Nh": 286.0,
    "Fl": 289.0,
    "Mc": 290.0,
    "Lv": 293.0,
    "Ts": 294.0,
    "Og": 2949.0,  # wrong?
}


def calculate_weight(chem: str) -> Dict[str, float]:
    # 匹配 Ti-22Al-23Nb-1Mo-1Zr 格式化学元素符号和相应的原子数
    ele_number_list = re.findall("([0-9]*)([A-Z][a-z]?)", chem)
    ele_number_list.reverse()
    # print(ele_number_list)

    weight_total = 0.0
    count_less = 0.0
    weight_list = []
    for atomic_percent, element in ele_number_list:
        if atomic_percent != "":
            atomic_percent = int(atomic_percent) * 0.01
            count_less += atomic_percent
        else:
            atomic_percent = 1 - count_less
        element_weight = atomic_mass_dict[element] * atomic_percent
        weight_list.append((element, element_weight))
        weight_total += element_weight

    ele_weight_dict = {
        element: round(percent / weight_total, 4) for element, percent in weight_list
    }
    print(ele_weight_dict)

    return ele_weight_dict


if __name__ == "__main__":
    chem = "Ti-22Al-23Nb-1Mo-1Zr"
    calculate_weight(chem)

# output:
# {'Zr': 0.0167, 'Mo': 0.0176, 'Nb': 0.3918, 'Al': 0.1088, 'Ti': 0.4651}
