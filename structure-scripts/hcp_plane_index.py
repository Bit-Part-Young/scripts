#!/usr/bin/env python3

"""
HCP 面指数的三指数和四指数坐标转换

示例:
(1, 1, 1) -> (1, 1, -2, 1)
(1, 0, 0) -> (1, 0, -1, 0)
(-1, -1, 0) -> (-1, -1, 2, 0)

(1, 1, -2, 0) -> (1, 1, 0)
(1, -1, 0, 0) -> (1, -1, 0)
(0, 1, -1, 0) -> (0, 1, 0)
"""

import argparse

import numpy as np


def f2t(four_index: list[int]) -> list[int]:
    """四指数坐标 -> 三指数坐标"""

    matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])

    three_index = matrix @ np.array(four_index).T

    return three_index.astype(int).tolist()


def t2f(three_index: list[int]) -> list[int]:
    """三指数坐标 -> 四指数坐标"""

    matrix = np.array([[1, 0, 0], [0, 1, 0], [-1, -1, 0], [0, 0, 1]])

    four_index = matrix @ np.array(three_index).T

    return four_index.astype(int).tolist()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Convert three index to four index")
    parser.add_argument(
        "-pi",
        "--plane_index",
        type=int,
        nargs="+",
        metavar="plane_index",
        help="plane index",
    )

    args = parser.parse_args()
    plane_index = args.plane_index

    if len(plane_index) == 3:
        four_index = t2f(plane_index)
        print(f"{tuple(plane_index)} -> {tuple(four_index)}")
    elif len(plane_index) == 4:
        three_index = f2t(plane_index)
        print(f"{tuple(plane_index)} -> {tuple(three_index)}")
