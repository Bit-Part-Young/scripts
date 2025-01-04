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

    matrix = np.array([[1, 0, -1, 0], [0, 1, -1, 0], [0, 0, 0, 1]])

    three_index = matrix @ np.array(four_index).T
    three_index = three_index.astype(int)

    if np.all(three_index % 3 == 0):
        three_index = three_index // 3
    elif np.all(three_index % 6 == 0):
        three_index = three_index // 6

    return three_index.tolist()


def t2f(three_index: list[int]) -> list[int]:
    """三指数坐标 -> 四指数坐标"""

    matrix = np.array([[2, -1, 0], [-1, 2, 0], [-1, -1, 0], [0, 0, 3]]) / 3

    four_index = matrix @ np.array(three_index).T

    # 非零
    four_index_tmp = four_index[np.abs(four_index) > 0.01]
    if np.any(np.abs(four_index_tmp) > 0.17) and np.any(np.abs(four_index_tmp) < 0.34):
        four_index = four_index * 3
    elif np.all(np.abs(four_index_tmp) < 0.17):
        four_index = four_index * 6

    four_index = four_index.astype(int)

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
