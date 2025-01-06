#!/usr/bin/env python3

"""
HCP 方向指数的三指数和四指数坐标转换

示例:
[1, 0, 0] -> [2, -1, -1, 0]
[0, 1, 0] -> [-1, 2, -1, 0]
[1, 1, 0] -> [1, 1, -2, 0]
[-1, -1, 0] -> [-1, -1, 2, 0]
[-1, 1, 0] -> [-1, 1, 0, 0]
[1, 0, 1] -> [2, -1, -1, 3]
[1, 1, 2] -> [1, 1, -2, 6]

[-1, -1, 2, 3] -> [-1, -1, 1]
[1, 0, -1, 0] -> [2, 1, 0]
"""

import argparse

import numpy as np


def f2t(four_index: list[int]) -> list[int]:
    """四指数坐标 -> 三指数坐标"""

    if len(four_index) != 4:
        raise ValueError("The length of four index must be 4. Exit!")
    elif sum(four_index[:-1]) != 0:
        raise ValueError("The sum of the first three index must be 0. Exit!")

    matrix = np.array([[1, 0, -1, 0], [0, 1, -1, 0], [0, 0, 0, 1]])

    three_index = matrix @ np.array(four_index).T

    if np.all(three_index % 3 == 0):
        three_index = three_index // 3

    return three_index.tolist()


def t2f(three_index: list[int]) -> list[int]:
    """三指数坐标 -> 四指数坐标"""

    if len(three_index) != 3:
        raise ValueError("The length of three index must be 3. Exit!")

    matrix = np.array([[2, -1, 0], [-1, 2, 0], [-1, -1, 0], [0, 0, 3]])

    four_index = matrix @ np.array(three_index).T

    if np.all(four_index % 3 == 0):
        four_index = four_index // 3

    return four_index.tolist()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Convert three index to four index")
    parser.add_argument(
        "-di",
        "--direction_index",
        type=int,
        nargs="+",
        metavar="direction_index",
        help="direction index",
    )

    args = parser.parse_args()
    direction_index = args.direction_index

    if len(direction_index) == 3:
        four_index = t2f(direction_index)
        print(f"{direction_index} -> {four_index}")
    elif len(direction_index) == 4:
        three_index = f2t(direction_index)
        print(f"{direction_index} -> {three_index}")
