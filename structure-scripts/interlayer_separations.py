#!/usr/bin/env python3

"""计算原子层间距变化（适用于表面/界面模型弛豫前后的原子层间距变化）"""

import argparse

import numpy as np
import pandas as pd
from ase.io import read


def get_interlayer_distance(
    structure_fn: str, precision: float = 0.001
) -> tuple[np.ndarray, int]:
    """统计原子层间距"""

    atoms = read(structure_fn)

    positions = atoms.positions

    z_coords = positions[:, 2]
    z_coords_rounded = precision * np.round(z_coords / precision)

    z_unique = np.unique(z_coords_rounded)

    distance = np.diff(z_unique)

    return distance, len(z_unique)


def interlayer_separations_cal(
    structure1_fn: str,
    structure2_fn: str,
    precision: float = 0.001,
) -> pd.DataFrame:
    """计算原子层间距变化（适用于表面/界面模型弛豫前后的原子层间距变化）"""

    distance1, layer1_count = get_interlayer_distance(structure1_fn, precision)
    distance2, layer2_count = get_interlayer_distance(structure2_fn, precision)

    if layer1_count != layer2_count:
        raise ValueError("Layer count of two structures is different! Exit.")
    else:
        interlayer_separations = distance1 - distance2
        ratio = (100 * (interlayer_separations / distance2)).round(2)

        print(f"Inerlayer separations info:")

        data = {
            "Interlayer_Distance_1": distance1,
            "Interlayer_Distance_2": distance2,
            "Interlayer_Separations": interlayer_separations,
            "Ratio(%)": ratio,
        }
        index = [f"{i}-{i+1}" for i in range(1, len(interlayer_separations) + 1)]

        interlayer_separations_info = pd.DataFrame(
            data=data,
            index=index,
        )

        # 逆序
        interlayer_separations_info = interlayer_separations_info[::-1]

    return interlayer_separations_info


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Interlayer separations calculation (for Surface/Interface model after and before relaxation).",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure1_fn",
        nargs="?",
        default="CONTCAR",
        metavar="STR",
        help="Structure filename after relaxation",
    )
    parser.add_argument(
        "structure2_fn",
        nargs="?",
        default="POSCAR",
        metavar="STR",
        help="Structure filename before relaxation",
    )

    args = parser.parse_args()

    info = interlayer_separations_cal(
        structure1_fn=args.structure1_fn, structure2_fn=args.structure2_fn
    )
    print(info)
