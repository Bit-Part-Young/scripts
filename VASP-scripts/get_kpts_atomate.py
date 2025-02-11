#!/usr/bin/env python3

"""获取 atomate 中不同 workflows 的默认 kpts 设置"""

import argparse

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints


def get_kpts_atomate(structure_fn: str = "POSCAR", force_gamma: bool = False):
    """获取 atomate 中不同 workflows 的默认 kpts 设置"""

    structure = Structure.from_file(structure_fn)

    print("Structure info:\n")
    print(structure)

    kpt_opt = Kpoints.automatic_density_by_vol(
        structure,
        kppvol=64,
        force_gamma=force_gamma,
    )

    kpt_static = Kpoints.automatic_density_by_vol(
        structure,
        kppvol=100,
        force_gamma=force_gamma,
    )

    kpt_elastic = Kpoints.automatic_density(
        structure,
        kppa=7000,
        force_gamma=force_gamma,
    )

    print("\nkpts info:\n")
    print(f"kpts for relax: {kpt_opt.kpts[0]}, style: {kpt_opt.style}")
    print(f"kpts for static: {kpt_static.kpts[0]}, style: {kpt_static.style}")
    print(f"kpts for elastic: {kpt_elastic.kpts[0]}, style: {kpt_elastic.style}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get defaults kpts of different workflows in atomate.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn",
        nargs="?",
        const="POSCAR",
        default="POSCAR",
        type=str,
        help="Structure filename",
    )

    parser.add_argument(
        "-fg",
        "--force_gamma",
        action="store_true",
        help="Force a gamma-centered mesh",
    )

    args = parser.parse_args()

    get_kpts_atomate(args.structure_fn)
