#!/usr/bin/env python3

"""获取 VASP 计算目录的数据"""

import sys

from atomate.vasp.drones import VaspDrone
from pymatgen.core.structure import Structure


def get_vasp_data(path: str) -> dict:
    """获取 VASP 计算目录的数据"""

    drone = VaspDrone()
    document = drone.assimilate(path=path)

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

    print(data_dict)

    return data_dict


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "."

    get_vasp_data(path)
