#!/usr/bin/env python3

"""获取 atomate VASP 计算目录的计算信息"""

import argparse
import pandas as pd
import os
from pymatgen.io.vasp import Vasprun, Outcar


def format_time(time_cost: float):
    hour = int(time_cost // 3600)
    minute = int((time_cost % 3600) // 60)
    second = int(time_cost % 60)

    return f"{hour:02d}h {minute:02d}m {second:02d}s"


def get_atomate_calc_info(root_dir: str, flag: bool = False):
    """获取 atomate VASP 计算目录的计算信息"""

    os.chdir(root_dir)

    # 仅当 添加 --flag 时，才显示保存的 csv 文件
    csv_fn = "atomate_calc_info.csv"
    if os.path.exists(csv_fn) and flag:
        if os.path.getsize(csv_fn) > 0.0:
            df = pd.read_csv(csv_fn)

            print()
            print(df)
    else:

        dir_list = []

        for dir_name in os.listdir("."):
            if os.path.isdir(dir_name):
                dir_list.append(dir_name)

        # 若计算完成，计算目录下的所有文件会被压缩
        outcar_fn_list = ["OUTCAR", "OUTCAR.gz", "OUTCAR.relax2.gz"]
        vasprun_fn_list = ["vasprun.xml", "vasprun.xml.gz", "vasprun.xml.relax2.gz"]

        data_list = []
        for dir_name in dir_list:
            # 首先找到该目录下的 FW-- 文件
            namefile = None
            for fn in os.listdir(dir_name):
                if "FW--" in fn:
                    namefile = fn
                    if "gz" in namefile:
                        namefile = namefile.replace(".gz", "")
                    break

            if namefile is None:
                continue

            # 然后尝试找到对应的 OUTCAR 和 vasprun.xml 文件
            found_files = False
            for outcar_fn, vasprun_fn in zip(outcar_fn_list, vasprun_fn_list):
                outcar_fn = os.path.join(dir_name, outcar_fn)
                vasprun_fn = os.path.join(dir_name, vasprun_fn)

                if os.path.exists(outcar_fn) and os.path.exists(vasprun_fn):
                    if "gz" in outcar_fn:
                        outcar = Outcar(outcar_fn)
                        time_cost = outcar.run_stats["Total CPU time used (sec)"]
                        time_cost = format_time(time_cost)

                        vasprun = Vasprun(vasprun_fn)
                        nsteps = int(vasprun.nionic_steps)
                        energy = vasprun.final_energy
                    else:
                        time_cost = None
                        nsteps = None
                        energy = None

                    data_dict = {
                        "Ion_Step": nsteps,
                        "Energy": energy,
                        "Time_Cost": time_cost,
                        "namefile": namefile,
                    }
                    data_list.append(data_dict)

                    found_files = True

                    break

            # 如果没有找到任何 OUTCAR/vasprun.xml 文件，添加一个只有 namefile 的记录
            if not found_files:
                data_dict = {
                    "Ion_Step": None,
                    "Energy": None,
                    "Time_Cost": None,
                    "namefile": namefile,
                }
                data_list.append(data_dict)

        df = pd.DataFrame(data_list)
        print()
        print(df)

        df.to_csv(csv_fn, index=False, sep=",")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root_dir", type=str, help="root directory")
    parser.add_argument("--flag", action="store_true", help="show saved csv file")

    args = parser.parse_args()

    get_atomate_calc_info(args.root_dir, args.flag)
