#!/usr/bin/env python3

"""将 VASP XDATCAR 文件转换为 NEP xyz 文件，并指定间隔抽样"""

import argparse

from ase.io import read, write


def xdatcar2xyz(
    xdatcar_fn: str = "XDATCAR",
    interval: int = 1,
    xyz_fn: str = "xdatcar.xyz",
):
    """将 VASP XDATCAR 文件转换为 NEP xyz 文件，并指定间隔"""

    atoms_list = read(
        xdatcar_fn,
        index=f"{interval - 1}::{interval}",
        format="vasp-xdatcar",
    )

    write(xyz_fn, atoms_list, format="extxyz", append=True)

    print(f"Total {len(atoms_list)} frames in {xdatcar_fn} converted to {xyz_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert VASP XDATCAR to NEP xzy file with specific interval sampling.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "xdatcar_fn",
        type=str,
        nargs="?",
        const="XDATCAR",
        default="XDATCAR",
        help="VASP XDATCAR filename.",
    )

    parser.add_argument(
        "xyz_fn",
        type=str,
        nargs="?",
        const="xdatcar.xyz",
        default="xdatcar.xyz",
        help="NEP xyz filename.",
    )

    parser.add_argument(
        "-i",
        "--interval",
        type=int,
        nargs="?",
        const=1,
        default=1,
        help="The interval of the frames to be converted.",
    )

    args = parser.parse_args()

    xdatcar2xyz(
        xdatcar_fn=args.xdatcar_fn,
        xyz_fn=args.xyz_fn,
        interval=args.interval,
    )
