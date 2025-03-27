#!/usr/bin/env python3

"""
构型文件格式转换

支持的输入格式: POSCAR, CONTCAR, XDATCAR, OUTCAR, vasprun.xml, *.vasp, *.xsd, *.xyz, *.lammps-data, ...

支持的输出格式: POSCAR, *.vasp, *.xsd, *.xyz, *.lammps-data, ...
"""

import argparse
import os

from ase.atoms import Atoms
from ase.io import read, write


def posconv(
    input_fn: str,
    output_fn: str,
):
    """构型文件格式转换"""

    input_fn_basename = os.path.basename(input_fn)
    input_format = input_fn_basename.split(".")[-1]

    # 解析 OUTCAR 有时会报错
    if input_format in ["XDATCAR", "xml", "OUTCAR"]:
        atoms: list[Atoms] = read(input_fn, index=":")
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
    elif output_format in ["lammps-data", "lmp"]:
        write_param_dict = {
            "format": "lammps-data",
            "atom_style": "atomic",
            "masses": True,
        }
    elif output_format in ["xyz", "extxyz"]:
        write_param_dict = {
            "format": "extxyz",
            "append": True,
        }

    write(
        output_fn,
        images=atoms,
        **write_param_dict,
    )

    print(f"\nConvert {input_fn} to {output_fn}!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Structure file format convert. Support most ASE recognized formats.",
        epilog="Author: SLY.",
        allow_abbrev=True,
    )

    parser.add_argument(
        "input_fn",
        type=str,
        help="Input filename, include: POSCAR, CONTCAR, XDATCAR、OUTCAR, vasprun.xml, *.vasp, *.xsd, *.xyz, *.lammps-data, ...",
    )

    parser.add_argument(
        "output_fn",
        type=str,
        help="Output filename, include: POSCAR, *.vasp, *.xsd, *.xyz, *.lammps-data, ...",
    )

    args = parser.parse_args()

    input_fn = args.input_fn
    output_fn = args.output_fn

    posconv(
        input_fn=input_fn,
        output_fn=output_fn,
    )
