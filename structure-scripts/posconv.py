#!/usr/bin/env python3

"""
构型文件格式转换

支持的输入格式: POSCAR, CONTCAR, XDATCAR, OUTCAR, vasprun.xml, *.vasp, *.xsd, *.xyz, *.lmp, *.lammps-data, ...

支持的输出格式: POSCAR, *.vasp, *.xsd, *.xyz, *.lmp, *.lammps-data, ...
"""

import argparse
import os

from ase.atoms import Atoms
from ase.build.tools import sort
from ase.io import read, write


def posconv(
    input_fn: str,
    output_fn: str,
    specorder: bool = False,
):
    """构型文件格式转换"""

    input_fn_basename = os.path.basename(input_fn)
    input_format = input_fn_basename.split(".")[-1]

    # 解析 OUTCAR 有时会报错
    if input_format in ["XDATCAR", "xml", "OUTCAR"]:
        atoms: list[Atoms] = read(input_fn, index=":")
        print(f"{input_fn} Frames: {len(atoms)}.")
    elif input_format in ["lmp", "lammps-data"]:
        atoms: Atoms = read(input_fn, format="lammps-data")
    elif input_format in ["xyz", "extxyz"]:
        atoms: list[Atoms] = read(input_fn, index=":", format="extxyz")
        # 对单个构型进行元素排序
        atoms = [sort(a) for a in atoms]
        print(f"{input_fn} Frames: {len(atoms)}.")
    else:
        atoms: Atoms = read(input_fn)

    output_fn_basename = os.path.basename(output_fn)
    output_format = output_fn_basename.split(".")[-1]
    if output_format in ["vasp", "POSCAR"]:
        write_param_dict = {
            "format": "vasp",
            "direct": True,
            "sort": True,
        }

    elif output_format in ["XDATCAR"]:
        write_param_dict = {
            "format": "vasp-xdatcar",
            "append": True,
        }

    elif output_format in ["lmp", "lammps-data"]:
        write_param_dict = {
            "format": "lammps-data",
            "atom_style": "atomic",
            "masses": True,
        }

        if specorder:
            # 按照给定元素顺序排序
            element_list = list(set(atoms.get_chemical_symbols()))
            element_sequence_list = ["Ti", "Al", "Nb", "Mo", "Zr", "V"]
            if all(element in element_sequence_list for element in element_list):
                write_param_dict["specorder"] = element_sequence_list

    elif output_format in ["xyz", "extxyz"]:
        write_param_dict = {
            "format": "extxyz",
            "append": True,
        }

    else:
        raise ValueError(f"Unsupported output format: {output_format}.")

    write(
        output_fn,
        images=atoms,
        **write_param_dict,
    )

    print(f"\nConvert {input_fn} to {output_fn}!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Structure file format convert. Support most ASE recognized formats.",
        allow_abbrev=True,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "input_fn",
        type=str,
        help="input filename, include: POSCAR, CONTCAR, XDATCAR, OUTCAR, vasprun.xml, *.vasp, *.xsd, *.xyz, *.lmp, *.lammps-data, ...",
    )

    parser.add_argument(
        "output_fn",
        type=str,
        help="output filename, include: POSCAR, *.vasp, *.xsd, *.xyz, *.lmp, *.lammps-data, ...",
    )

    parser.add_argument(
        "--specorder",
        action="store_true",
        help="whether to enable speorder parameter in LAMMPS data format.",
    )

    args = parser.parse_args()

    posconv(
        input_fn=args.input_fn,
        output_fn=args.output_fn,
        specorder=args.specorder,
    )
