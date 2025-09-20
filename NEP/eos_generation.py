#!/usr/bin/env python3

"""生成 BCC/FCC/HCP 的 EOS 构型（晶格常数附近 0.2 Å 范围内，体积应变范围通常宽于 -10%~10%）"""

import os
import argparse

import numpy as np
from ase.build import bulk
from ase.io import write


def eos_configurations_cubic_generation(element: str, a: float, crystalstructure: str):

    structure_folder = "0-structures"
    os.makedirs(structure_folder, exist_ok=True)

    a = round(a, 2)
    flag = 0
    for i, lc in enumerate(np.arange(a - 0.2, a + 0.2, 0.01), start=1):
        lc = round(lc, 4)
        if crystalstructure in ["fcc", "bcc"]:
            atoms = bulk(element, crystalstructure, a=lc, cubic=True)

        structure_fn = f"{structure_folder}/POSCAR.{i}"
        write(structure_fn, atoms, format="vasp", direct=True, sort=True, vasp5=True)
        flag += 1

    print(f"Total {flag} Configurations generated!")


def eos_configurations_hcp_generation(
    element: str, a: float, c: float, crystalstructure: str = "hcp"
):
    covera = round(c / a, 3)
    a_initial = a
    a = round(a, 2)

    # --------------------- a 变化 c/a 不变 ---------------------
    folder = f"0-structures-a"
    os.makedirs(folder, exist_ok=True)

    scale = 0.2
    flag = 0
    for i, lc in enumerate(
        np.arange(round(a - scale, 2), round(a + scale + 0.01, 2), 0.01), start=1
    ):
        lc = round(lc, 4)
        atoms = bulk(element, crystalstructure, a=lc, covera=covera)

        structure_fn = f"{folder}/POSCAR.{i}"
        write(structure_fn, atoms, format="vasp", direct=True, sort=True, vasp5=True)

        flag += 1

    print(f"Total {flag} Configurations with a variation generated!")

    # --------------------- c/a 变化 a 不变 ---------------------

    folder = f"0-structures-ca"
    os.makedirs(folder, exist_ok=True)

    a = a_initial
    scale = 0.1
    covera = round(covera, 2)

    flag = 0
    for i, lc in enumerate(
        np.arange(round(covera - scale, 2), round(covera + scale + 0.01, 2), 0.01),
        start=1,
    ):
        lc = round(lc, 4)
        atoms = bulk(element, crystalstructure, a=a, covera=lc)

        structure_fn = f"{folder}/POSCAR.{i}"
        write(structure_fn, atoms, format="vasp", direct=True, sort=True, vasp5=True)

        flag += 1

    print(f"Total {flag} Configurations with c/a variation generated!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate EOS configurations of BCC/FCC/HCP with lattice constants variation in -0.2~0.2 Å.",
        epilog="Author: SLY.",
    )
    parser.add_argument("element", type=str, help="element symbol")
    parser.add_argument(
        "crystalstructure",
        type=str,
        choices=["bcc", "fcc", "hcp"],
        help="crystal structure",
    )
    parser.add_argument(
        "-lc",
        "--lattice_constants",
        type=float,
        nargs="+",
        metavar="FLOAT",
        help="lattice constants",
    )
    args = parser.parse_args()

    element = args.element
    crystalstructure = args.crystalstructure
    lattice_constants = args.lattice_constants
    if len(lattice_constants) == 1:
        a = lattice_constants[0]
        c = a
    elif len(lattice_constants) == 2:
        a, c = lattice_constants
    else:
        raise ValueError(f"-lc, --lattice_constants allows only 1 or 2 arguments.")

    if crystalstructure in ["fcc", "bcc"]:
        eos_configurations_cubic_generation(element, a, crystalstructure)
    elif crystalstructure == "hcp":
        eos_configurations_hcp_generation(
            element, a=a, c=c, crystalstructure=crystalstructure
        )
    else:
        raise ValueError(f"Invalid crystal structure: {crystalstructure}.")
