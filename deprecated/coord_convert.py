#!/usr/bin/env python3

"""
基于 pymatgen 进行 Cartesian 和 Direct 坐标转换
Deprecated
"""

import os
import sys
from fnmatch import fnmatch

from print_help import print_help
from pymatgen.io.vasp.inputs import Poscar


def get_frac_flag(poscar_fn="POSCAR"):
    """
    get coordinates format conversion flag from POSCAR file.
    """

    with open(poscar_fn, "r") as f:
        lines = f.readlines()
        cha = lines[7].strip()
        comment = lines[0].strip()

        frac_flag = True
        if cha.startswith("d") or cha.startswith("D"):
            frac_flag = False

        return frac_flag, comment


def write_poscar(poscar_fn="POSCAR"):
    """
    write to POSCAR file.
    """

    frac_flag, comment = get_frac_flag(poscar_fn=poscar_fn)

    poscar = Poscar.from_file(poscar_fn)
    poscar.write_file(poscar_fn, direct=frac_flag)

    os.popen(f"sed -i '1s/.*/{comment}/g' {poscar_fn}")

    if frac_flag:
        print(f"The {poscar_fn} file is converted to direct coordinates.")
    else:
        print(f"The {poscar_fn} file is converted to cartesian coordinates.")


def main():
    argv_dict = {
        "poscar_file": "POSCAR file, supported format: *POSCAR*, *.vasp, *.poscar"
    }
    description = "POSCAR file coordinates format conversion calling pymatgen."

    poscar_fn = sys.argv[1]

    if sys.argv[1] in ["-h", "--help"]:
        print_help(argv_dict, description)
    elif len(sys.argv) > 2:
        print(f"\nError: except 1 argument, but get {len(sys.argv)-1} arguments.")
        print_help(argv_dict, description)
    elif len(sys.argv) == 2 and (
        fnmatch(poscar_fn, "*POSCAR*")
        or fnmatch(poscar_fn, "*.vasp")
        or fnmatch(poscar_fn, "*.poscar")
    ):
        write_poscar(poscar_fn=poscar_fn)


if __name__ == "__main__":
    main()
