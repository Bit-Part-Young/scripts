#!/usr/bin/env python3

"""
根据 material_id 检查、获取 MP phonon 数据并绘制band structure
经简要测试，只有 mp-149 (Si) 有 MP phonon 数据，mp-46 mp-129 等没有
"""

import argparse
import os

from mp_api.client import MPRester
from pymatgen.phonon.plotter import PhononBSPlotter

API_KEY = os.getenv("PMG_MAPI_KEY")
if API_KEY is None:
    raise OSError("PMG_MAPI_KEY environment variable is not set. Please check!")


def check(material_id: str) -> bool:
    with MPRester(API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=material_id, has_props=["phonon"]
        )

    return len(docs) > 0


def get_phonon_data(material_id):

    with MPRester(API_KEY) as mpr:
        phonon_doc = mpr.materials.phonon.search(material_ids=material_id)

    return phonon_doc


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check, get and plot MP phonon data.", epilog="Author: SLY."
    )

    parser.add_argument("material_id", help="MP material id")

    args = parser.parse_args()

    if not check(args.material_id):
        print(f"Material {args.material_id} not found")
        exit(1)

    phonon_doc = get_phonon_data(args.material_id)

    phonon_data = phonon_doc[0]
    ph_bs = phonon_data.ph_bs
    ph_dos = phonon_data.ph_dos

    phonon_bs_plotter = PhononBSPlotter(bs=ph_bs)

    fig_fn = f"phonon_bs_{args.material_id}.png"
    phonon_bs_plotter.save_plot(fig_fn)

    print(f"\nPhonon band structure plot saved to {fig_fn}.")
