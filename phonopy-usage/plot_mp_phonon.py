#!/usr/bin/env python3

"""
根据 material_id 获取 Materials Project 声子数据并绘制 band structure

API key 取自环境变量 PMG_MAPI_KEY。MP 的声子方法有 dfpt、pheasy 两种，部分材料两者都有；
未指定 --method 时用第一个可用的。没有声子数据的材料（如 mp-75 Nb）会报错退出。
"""

import argparse
import sys

from mp_api.client import MPRester
from mp_api.client.core.exceptions import MPRestError
from pymatgen.phonon.plotter import PhononBSPlotter

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get and plot MP phonon band structure.", epilog="Author: SLY."
    )
    parser.add_argument("material_id", help="MP material id, e.g. mp-149")
    parser.add_argument("-m", "--method", help="phonon method: dfpt or pheasy (default: first available)")
    parser.add_argument("-o", "--output_fn", help="output figure (default: phonon_bs_<id>_<method>.png)")
    args = parser.parse_args()

    with MPRester(mute_progress_bars=True) as mpr:
        try:
            docs = mpr.materials.phonon.search(material_ids=[args.material_id], fields=["phonon_method"])
        except MPRestError as e:
            sys.exit(str(e))
        methods = [d.phonon_method.value for d in docs]
        print(f"{args.material_id}: phonon methods {methods}")

        method = args.method or methods[0]
        if method not in methods:
            sys.exit(f"no {method} phonon data for {args.material_id}")
        bs = mpr.materials.phonon.get_bandstructure_from_material_id(args.material_id, phonon_method=method)

    fig_fn = args.output_fn or f"phonon_bs_{args.material_id}_{method}.png"
    PhononBSPlotter(bs=bs.to_pmg).save_plot(fig_fn)

    print(f"Phonon band structure plot saved to {fig_fn}")
