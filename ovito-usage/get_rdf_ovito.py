#!/usr/bin/env python3

"""
计算轨迹文件平均 RDF
reference: https://gitee.com/mayuan_JLUPHY/my_script/blob/master/lammps-toolkit/getRDFovito.py
"""

import argparse

from ovito.io import export_file, import_file
from ovito.modifiers import CoordinationAnalysisModifier, TimeAveragingModifier
from ovito.pipeline import Pipeline


def rdf_cal(
    trajectory_fn: str,
    cutoff: float,
    num_bins: int,
    output_fn: str,
    partial: bool,
):

    pipeline = import_file(trajectory_fn)
    pipeline: Pipeline
    print("Number of MD frames:", pipeline.num_frames)

    # 输出体系的原子种类及其 ID
    for type in pipeline.compute().particles.particle_types.types:
        print(f"Type {type.id}: {type.name}")

    # 添加  RDF 计算 modifier
    pipeline.modifiers.append(
        CoordinationAnalysisModifier(
            cutoff=cutoff,
            number_of_bins=num_bins,
            partial=partial,
        )
    )

    # 对所有帧的 RDF 数据进行时间平均
    pipeline.modifiers.append(TimeAveragingModifier(operate_on="table:coordination-rdf"))

    # 导出
    export_file(
        data=pipeline,
        file=output_fn,
        format="txt/table",
        key="coordination-rdf[average]",
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate (average) RDF with ovito Python package.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
    )

    parser.add_argument(
        "-i",
        "--input_fn",
        default="XDATCAR",
        required=True,
        metavar="trajectory_fn",
        help="configuration/trajectory filename",
    )

    parser.add_argument(
        "-c",
        "--cutoff",
        default=5.0,
        type=float,
        metavar="cutoff",
        help="cutoff radius",
    )

    parser.add_argument(
        "-n",
        "--num_bins",
        default=100,
        type=int,
        metavar="num_bins",
        help="number of bins",
    )

    parser.add_argument(
        "-o",
        "--output_fn",
        default="rdf.txt",
        type=str,
        metavar="output_fn",
        help="output filename",
    )

    parser.add_argument(
        "-p",
        "--partial",
        action="store_true",
        help="whether to get partial element RDF",
    )

    args = parser.parse_args()

    input_fn = args.input_fn
    cutoff = args.cutoff
    num_bins = args.num_bins
    output_fn = args.output_fn
    partial = args.partial

    rdf_cal(
        trajectory_fn=input_fn,
        cutoff=cutoff,
        num_bins=num_bins,
        output_fn=output_fn,
        partial=partial,
    )
