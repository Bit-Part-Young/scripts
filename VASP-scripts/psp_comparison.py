#!/usr/bin/env python3

"""比较 VASP 和 pymatgen 推荐的常用元素 PBE 赝势 (VASP.5.4.4)"""

import argparse
import os
import subprocess

import pandas as pd
from monty.os.path import zpath
from pymatgen.core import SETTINGS

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
}


def get_shell_outputs(command: str) -> str:
    """获取 Shell 命令的输出"""

    result = subprocess.run(
        command,
        shell=True,
        text=True,
        capture_output=True,
    )
    return result.stdout


def get_potcar_info(psp_symbol: str):
    """获取 POTCAR 文件信息"""

    PMG_VASP_PSP_DIR = SETTINGS.get("PMG_VASP_PSP_DIR")

    functional_subdir = "POT_GGA_PAW_PBE"

    paths_to_try: list[str] = [
        os.path.join(PMG_VASP_PSP_DIR, functional_subdir, f"POTCAR.{psp_symbol}"),
        os.path.join(PMG_VASP_PSP_DIR, functional_subdir, psp_symbol, "POTCAR"),
    ]

    for path in paths_to_try:
        path = os.path.expanduser(path)
        path = zpath(path)
        if os.path.isfile(path):

            command = f"head -n 1 {path}"
            titel = get_shell_outputs(command).strip()

            command = f"grep ENMAX {path} | awk '{{print $3}}'"
            enmax = float(get_shell_outputs(command).strip()[:-1])

            command = f"grep ZVAL {path} | awk '{{print $6}}'"
            zval = int(float(get_shell_outputs(command).strip()))

            command = f"grep VRHFIN {path} | awk -F':' '{{print $2}}'"
            vrhfin = get_shell_outputs(command).strip()

            psp_info_dict = {
                "TITEL": titel,
                "ENMAX": enmax,
                "ZVAL": zval,
                "VRHFIN": vrhfin,
            }

    return psp_info_dict


def potcar_pbe_compare(element_symbol: str):
    """比较 VASP 和 pymatgen 推荐的常用元素 PBE 赝势 (VASP.5.4.4)"""

    psp_vasp = pbe_dict_vasp[element_symbol]
    psp_pymatgen = pbe_dict_pymatgen[element_symbol]

    psp_info_dict_vasp = get_potcar_info(psp_symbol=psp_vasp)
    psp_info_dict_pymatgen = get_potcar_info(psp_symbol=psp_pymatgen)
    if psp_vasp != psp_pymatgen:
        print("VASP and pymatgen recommended is different.\n")
        df = pd.DataFrame(
            [psp_info_dict_vasp, psp_info_dict_pymatgen],
            index=["VASP", "pymatgen"],
        )
    else:
        print("VASP and pymatgen recommended is same.\n")
        df = pd.DataFrame(
            [psp_info_dict_vasp],
            index=["VASP/pymatgen"],
        )

    # 索引行输出宽度设置为 15
    df.index = df.index.map(lambda x: f"{x:<15}")

    print(df)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Compare VASP, pymatgen recommended PBE pseudopotentials.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "element_symbol",
        help="Element symbol, eg. Ti",
    )

    args = parser.parse_args()

    potcar_pbe_compare(element_symbol=args.element_symbol)
