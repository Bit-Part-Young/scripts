#!/usr/bin/env python3

"""获取轨迹文件中每帧构型的点缺陷（空位、间隙）数目"""

import argparse

from ovito.io import export_file, import_file
from ovito.modifiers import WignerSeitzAnalysisModifier
from ovito.pipeline import Pipeline


def get_point_defects(
    trajectory_fn: str,
    output_fn: str = "pd_count.txt",
):
    """获取轨迹文件中每帧构型的点缺陷（空位、间隙）数目"""

    pipeline: Pipeline = import_file(trajectory_fn)

    nframes = pipeline.num_frames

    if nframes < 2:
        raise ValueError("The trajectory file should contain at least 2 frames.")
    else:
        print("Number of MD frames:", pipeline.num_frames)

    # 输出体系的原子种类及其 ID
    # for type in pipeline.compute().particles.particle_types.types:
    #     print(f"Type {type.id}: {type.name}")

    wsa_modifier = WignerSeitzAnalysisModifier()

    pipeline.modifiers.append(wsa_modifier)

    export_file(
        pipeline,
        output_fn,
        "txt/attr",
        multiple_frames=True,
        columns=["Frame", "WignerSeitz.vacancy_count", "WignerSeitz.interstitial_count"],
    )

    print(f"\nData is save to {output_fn}.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get number of point defects (vacancy & interstitial) from a trajectory file.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "trajectory_fn",
        type=str,
        metavar="trajectory_fn",
        default="dump.lammpstrj",
        help="The trajectory filename",
    )

    parser.add_argument(
        "-o",
        "--output_fn",
        type=str,
        metavar="output_fn",
        default="pd_count.txt",
        help="output data filename",
    )

    args = parser.parse_args()

    trajectory_fn = args.trajectory_fn
    output_fn = args.output_fn

    get_point_defects(
        trajectory_fn=trajectory_fn,
        output_fn=output_fn,
    )
