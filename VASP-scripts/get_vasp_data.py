#!/usr/bin/env python3

"""获取 VASP 计算目录的数据"""

import argparse

from atomate.vasp.drones import VaspDrone
from pymatgen.core.structure import Structure


def get_vasp_data(vasp_folder: str) -> dict:
    """获取 VASP 计算目vasp_folder"""

    drone = VaspDrone()
    document = drone.assimilate(path=vasp_folder)

    energy = document["output"]["energy"]

    structure_dict = document["output"]["structure"]
    structure = Structure.from_dict(structure_dict)
    lattice = structure.lattice.abc

    data_dict = {
        "energy": energy,
        "a": lattice[0],
        "b": lattice[1],
        "c": lattice[2],
    }

    return data_dict


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get data from VASP calculation directory with atomate package.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "vasp_folder",
        nargs="?",
        type=str,
        default=".",
        help="VASP calculation directory.",
    )

    args = parser.parse_args()

    vasp_folder = args.path

    data_dict = get_vasp_data(vasp_folder)

    print(data_dict)
