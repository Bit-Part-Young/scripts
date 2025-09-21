#!/usr/bin/env python3

"""将 VASP 输出文件 OUTCAR/vasprun.xml 转换为 xyz 文件"""

import argparse

from ase.io import read, write


def vaspout2xyz(
    vaspout_fn: str = "OUTCAR",
    xyz_fn: str = "vaspout.xyz",
    multiple_frame: bool = False,
):
    """将 VASP 输出文件 OUTCAR/vasprun.xml 转换为 xyz 文件"""

    if multiple_frame:
        atoms_list = read(vaspout_fn, index=":")
    else:
        atoms_list = [read(vaspout_fn, index=-1)]

    for atoms in atoms_list:
        virial = -atoms.get_volume() * atoms.get_stress(voigt=False)
        atoms.info["virial"] = virial.round(5)

        write(xyz_fn, atoms, format="extxyz", append=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert VASP outputs to xyz.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "vaspout_fn", type=str, help="VASP output filename, eg. OUTCAR / vasprun.xml"
    )
    parser.add_argument("xyz_fn", type=str, help="xyz filename")
    parser.add_argument(
        "-mf",
        "--multiple_frame",
        action="store_true",
        help="whether get multiple frames",
    )

    args = parser.parse_args()

    vaspout2xyz(
        vaspout_fn=args.vaspout_fn,
        xyz_fn=args.xyz_fn,
        multiple_frame=args.multiple_frame,
    )
