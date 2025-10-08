#!/usr/bin/env python3

"""
统计轨迹文件中的结构类型
可用于检查 SLC+NPH 中最后弛豫过程中轨迹是否全为液相
"""

import argparse
import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

import ovito._extensions.pyscript
import pandas as pd
from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier
from ovito.pipeline import Pipeline


def count_structure_types(trajectory_fn: str):
    """统计轨迹文件中的结构类型"""

    pipeline = import_file(trajectory_fn)
    pipeline: Pipeline
    nframes = pipeline.num_frames

    print("\nNumber of MD frames:", nframes)
    print()

    cna_modifier = CommonNeighborAnalysisModifier()
    pipeline.modifiers.append(cna_modifier)

    structure_type_keys = [
        "CommonNeighborAnalysis.counts.OTHER",
        "CommonNeighborAnalysis.counts.BCC",
        "CommonNeighborAnalysis.counts.FCC",
        "CommonNeighborAnalysis.counts.HCP",
        "CommonNeighborAnalysis.counts.ICO",
    ]

    count_list = []
    for frame in range(nframes):
        data = pipeline.compute(frame)
        attributes = data.attributes

        count_dict = {key: attributes.get(key, None) for key in structure_type_keys}
        count_list.append(count_dict)

    df = pd.DataFrame(count_list)
    df.columns = ["OTHER", "BCC", "FCC", "HCP", "ICO"]

    print(df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count the structure types in the trajectory with OVITO.",
        epilog="Author: SLY.",
    )

    parser.add_argument("trajectory_fn", help="trajectory filename")

    args = parser.parse_args()

    count_structure_types(args.trajectory_fn)
