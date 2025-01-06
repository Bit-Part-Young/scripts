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

    if len(four_index) != 4:
        raise ValueError("The length of four index must be 4. Exit!")
    elif sum(four_index[:-1]) != 0:
        raise ValueError("The sum of the first three index must be 0. Exit!")

    matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])

    three_index = matrix @ np.array(four_index).T

    return three_index.tolist()


def t2f(three_index: list[int]) -> list[int]:
    """三指数坐标 -> 四指数坐标"""

    if len(three_index) != 3:
        raise ValueError("The length of three index must be 3. Exit!")

    matrix = np.array([[1, 0, 0], [0, 1, 0], [-1, -1, 0], [0, 0, 1]])

    four_index = matrix @ np.array(three_index).T

    return four_index.tolist()


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
