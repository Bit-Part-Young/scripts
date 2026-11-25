#!/usr/bin/env python3

"""获取 VASP 计算目录输出数据（优先用 atomate VaspDrone，未安装时退回 pymatgen）

atomate、pymatgen 均为可选依赖、按需导入；都没有时请用 get_vasp_data2.py（仅需 pandas）。
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def seconds_to_hms(seconds: float):
    """将秒数转换为时分秒格式"""

    hours, remainder = divmod(seconds, 3600)
    minutes, seconds_new = divmod(remainder, 60)

    return f"{int(hours)}h-{int(minutes)}min-{int(seconds_new)}s"


def _read_atomate(path):
    from atomate.vasp.drones import VaspDrone
    from pymatgen.core.structure import Structure

    document = VaspDrone().assimilate(path=path)
    structure = Structure.from_dict(document["output"]["structure"])
    elapsed = document["run_stats"]["overall"]["Elapsed time (sec)"]
    return document["output"]["energy"], structure, elapsed


def _read_pymatgen(path):
    from pymatgen.io.vasp import Outcar, Vasprun

    vasprun = Vasprun(Path(path) / "vasprun.xml", parse_dos=False, parse_eigen=False)
    elapsed = Outcar(Path(path) / "OUTCAR").run_stats["Elapsed time (sec)"]
    return vasprun.final_energy, vasprun.final_structure, elapsed


def get_vasp_data(path: str = ".") -> pd.DataFrame:
    """获取 VASP 计算目录输出数据"""

    try:
        import atomate  # noqa: F401

        reader = _read_atomate
    except ImportError:
        try:
            import pymatgen  # noqa: F401
        except ImportError:
            sys.exit("Error: neither atomate nor pymatgen is installed; use get_vasp_data2.py instead.")
        reader = _read_pymatgen

    energy, structure, time_cost = reader(path)
    natoms = len(structure)
    time_cost = seconds_to_hms(time_cost)

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
        description="Get data from VASP calculation path with atomate (or pymatgen if atomate is missing).",
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
