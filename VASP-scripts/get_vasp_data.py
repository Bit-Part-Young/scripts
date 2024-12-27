#!/usr/bin/env python3

"""获取 VASP 计算目录输出数据（利用 atomate package）"""

import argparse

import pandas as pd
from atomate.vasp.drones import VaspDrone
from pymatgen.core.structure import Structure


def seconds_to_hms(seconds: float):
    """将秒数转换为时分秒格式"""

    hours, remainder = divmod(seconds, 3600)
    minutes, seconds_new = divmod(remainder, 60)

    return f"{int(hours)}h-{int(minutes)}min-{int(seconds_new)}s"


def get_vasp_data(path: str = ".") -> dict:
    """获取 VASP 计算目录输出数据"""

    drone = VaspDrone()
    document = drone.assimilate(path=path)

    energy = document["output"]["energy"]
    natoms = document["nsites"]
    time_cost = document["run_stats"]["overall"]["Total CPU time used (sec)"]
    time_cost = seconds_to_hms(time_cost)

    structure_dict = document["output"]["structure"]
    structure = Structure.from_dict(structure_dict)
    lattice = structure.lattice.abc

    data_dict = {
        "a": lattice[0],
        "b": lattice[1],
        "c": lattice[2],
        "natoms": natoms,
        "energy": energy,
        "time_cost": time_cost,
    }

    df = pd.DataFrame(data_dict, index=[0]).round(5)

    return df


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get data from VASP calculation path with atomate package.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "path",
        nargs="?",
        type=str,
        default=".",
        help="VASP calculation path",
    )

    args = parser.parse_args()

    path = args.path

    df = get_vasp_data(path)

    print(df)
