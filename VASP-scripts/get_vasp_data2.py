#!/usr/bin/env python3

"""获取 VASP 计算目录输出数据（利用 pymatgen package）"""

import argparse
import os

import pandas as pd
from pymatgen.io.vasp.outputs import Outcar, Vasprun


def seconds_to_hms(seconds: float):
    """将秒数转换为时分秒格式"""

    hours, remainder = divmod(seconds, 3600)
    minutes, seconds_new = divmod(remainder, 60)

    return f"{int(hours)}h-{int(minutes)}min-{int(seconds_new)}s"


def get_vasp_data(path: str = "."):
    """获取 VASP 计算目录输出数据"""

    vasprun_path = os.path.join(path, "vasprun.xml")
    outcar_path = os.path.join(path, "OUTCAR")

    vasprun = Vasprun(vasprun_path)
    outcar = Outcar(outcar_path)

    structure = vasprun.final_structure
    natoms = structure.num_sites

    time_cost = outcar.run_stats["Total CPU time used (sec)"]
    time_cost = seconds_to_hms(time_cost)

    data_dict = {
        "energy": vasprun.final_energy,
        "natoms": natoms,
        "converged": vasprun.converged,
        "time_cost": time_cost,
    }

    df = pd.DataFrame(data_dict, index=[0]).round(5)

    return df


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get data from VASP calculation path with pymatgen package.",
        epilog="Author: SLY",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "path",
        type=str,
        default=".",
        nargs="?",
        help="VASP calculation path",
    )

    args = parser.parse_args()

    path = args.path

    df = get_vasp_data(path)

    print(df)
