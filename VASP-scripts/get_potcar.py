#!/usr/bin/env python3

"""生成 VASP、pymatgen 推荐的赝势文件"""

import os
import sys

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Potcar

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


def get_potcar_pymatgen(
    structure_fn: str = "POSCAR",
    recommended: str = "vasp",
):
    """生成 VASP、pymatgen 推荐的赝势文件"""

    if os.path.exists("POTCAR"):
        os.remove("POTCAR")

    structure = Structure.from_file(structure_fn)

    element_list = structure.composition.as_dict().keys()

    if recommended == "vasp":
        psp_list = [vasp_pbe_dict[element] for element in element_list]

        print("\nPOTCAR recommended by VASP generated.")
    elif recommended == "pymatgen":
        psp_list = [pymatgen_pbe_dict[element] for element in element_list]

        print("\nPOTCAR recommended by pymatgen generated.")

    # 这里的 symbols 是元素赝势符号，而非元素符号
    potcar = Potcar(
        symbols=psp_list,
        functional="PBE",
    )
    potcar.write_file("POTCAR")


if __name__ == "__main__":

    structure_fn = sys.argv[1] if len(sys.argv) > 1 else "POSCAR"
    recommended = sys.argv[2] if len(sys.argv) > 2 else "vasp"

    get_potcar_pymatgen(
        structure_fn=structure_fn,
        recommended=recommended,
    )
