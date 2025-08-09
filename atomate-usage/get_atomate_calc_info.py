#!/usr/bin/env python3

"""获取 atomate VASP 计算目录的计算信息"""

import argparse
import pandas as pd
import os
from pymatgen.io.vasp import Vasprun, Outcar
import subprocess


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
        try:
            df = pd.read_csv(csv_fn)
            print()
            print(df)
        except pd.errors.EmptyDataError:
            print(f"\n{csv_fn} is empty.")
            exit()
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
            namefile = None
            for fn in os.listdir(dir_name):
                if "FW--" in fn:
                    namefile = fn
                    if "gz" in namefile:
                        namefile = namefile.replace(".gz", "")
                    break

            if namefile is None:
                continue

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
                        state = "Completed"
                    else:
                        try:
                            results = subprocess.run(
                                f"vasp_timecost.sh {dir_name}",
                                shell=True,
                                capture_output=True,
                            )

                            results_str = results.stdout.decode("utf-8").strip()

                            nsteps = int(results_str.split(" ")[0])
                            time_cost = results_str.split("cost: ")[1]
                            energy = round(
                                float(results_str.split("energy: ")[1].split(" ")[0]), 5
                            )
                        except:
                            print(f"vasp_timecost.sh excutable script not found!")

                            nsteps = None
                            time_cost = None
                            energy = None

                        state = "Running"

                    data_dict = {
                        "Ion_Step": nsteps,
                        "Energy": energy,
                        "State": state,
                        "Time_Cost": time_cost,
                        "namefile": namefile,
                        "launch_dir": dir_name,
                    }
                    data_list.append(data_dict)

                    found_files = True

                    break

            if not found_files:
                data_dict = {
                    "Ion_Step": None,
                    "Energy": None,
                    "Time_Cost": None,
                    "namefile": namefile,
                    "launch_dir": dir_name,
                    "State": None,
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
