#!/usr/bin/env python3

"""生成 VASP、pymatgen 推荐的 PBE 赝势 POTCAR 文件（不依赖 pymatgen，速度更快）"""

import argparse
import os
import subprocess

# Data source: https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/io/vasp/MPRelaxSet.yaml
pbe_dict_pymatgen = {
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
    "P": "P",
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
    "Pt": "Pt",
    "Au": "Au",
    "Pb": "Pb_d",
    "Ag": "Ag",
    "Pd": "Pd",
}

# Data source: https://www.vasp.at/wiki/index.php/Available_PAW_potentials
pbe_dict_vasp = {
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
    "P": "P",
    "S": "S",
    "Cl": "Cl",
    "Ca": "Ca_sv",
    "Sc": "Sc_sv",
    "Ti": "Ti_sv",
    "V": "V_sv",
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
    "Pt": "Pt",
    "Au": "Au",
    "Pb": "Pb_d",
    "Ag": "Ag",
    "Pd": "Pd",
}


def parse_poscar_file(structure_fn: str = "POSCAR") -> list:
    """解析 poscar 文件"""

    with open(structure_fn, "r") as file:
        lines = file.readlines()

    element_symbols = lines[5].strip().split()

    return element_symbols


def get_psp(
    structure_fn: str = "POSCAR",
    psp_recommended: str = "vasp",
):
    """生成 VASP、pymatgen 推荐的 PBE 赝势 POTCAR 文件"""

    if os.path.exists("POTCAR"):
        os.remove("POTCAR")

    psp_path_root = os.getenv("VASP_PSP_PATH")
    functional_subdir = "PAW_PBE"
    psp_path = os.path.join(psp_path_root, functional_subdir)

    element_symbols = parse_poscar_file(structure_fn)

    if psp_recommended == "vasp":
        psp_dict = pbe_dict_vasp

        print("\nPOTCAR recommended by VASP generated.")
    elif psp_recommended == "pymatgen":
        psp_dict = pbe_dict_pymatgen

        print("\nPOTCAR recommended by pymatgen generated.")

    psp_symbols = [psp_dict[element] for element in element_symbols]

    for psp_symbol in psp_symbols:
        element_psp_path = os.path.join(psp_path, psp_symbol)
        element_psp_fn = os.path.join(element_psp_path, "POTCAR")

        subprocess.run(f"cat {element_psp_fn} >> POTCAR", shell=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate VASP, pymatgen recommended PBE pseudopotentials (pymatgen independent, faster).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
        allow_abbrev=True,
    )

    parser.add_argument(
        "-pr",
        "--psp_recommended",
        nargs="?",
        const="vasp",
        default="vasp",
        type=str,
        choices=["vasp", "pymatgen"],
        help="recommended pseudopotential type",
    )

    parser.add_argument(
        "structure_fn",
        nargs="?",
        default="POSCAR",
        type=str,
        help="structure filename",
    )

    args = parser.parse_args()

    structure_fn = args.structure_fn
    psp_recommended = args.psp_recommended

    get_psp(
        structure_fn=structure_fn,
        psp_recommended=psp_recommended,
    )
