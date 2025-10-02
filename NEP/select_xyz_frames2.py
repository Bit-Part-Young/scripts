#!/usr/bin/env python3

"""
按比例选取指定个数的小于受力阈值的构型（用于缺陷、层错构型弛豫轨迹构型的筛选）
需借助 plot_energy_force.sh 脚本生成 force_info.dat 文件
"""

import argparse
import os

import numpy as np


def parse_xyz_file(input_xyz_fn: str) -> list:
    """解析 extxyz 文件"""

    with open(input_xyz_fn, "r") as file:
        lines = file.readlines()

    frames = []
    i = 0
    while i < len(lines):
        natoms = int(lines[i].strip())
        frame_info = lines[i + 1].strip()

        atoms_info = lines[i + 2 : i + 2 + natoms]
        frames.append((natoms, frame_info, atoms_info))
        i += 2 + natoms

    return frames


def write_xyz_file(frames: list, output_xyz_fn: str):
    """写入 extxyz 文件"""

    with open(output_xyz_fn, "w") as file:
        for natoms, frame_info, atoms_info in frames:
            file.write(f"{natoms}\n")
            file.write(f"{frame_info}\n")
            for atom_info in atoms_info:
                file.write(f"{atom_info.strip()}\n")


def force_exceeds_threshold(
    force_threshold: float = 5.0, num_selected: int = 5
) -> list:
    """按比例选取指定个数的小于受力阈值的构型索引"""

    if not os.path.exists("force_info.dat"):
        print(f"force_info.dat does not exist! Please check.")
        exit()

    force = np.loadtxt("force_info.dat", skiprows=1)

    if force.size == 1:
        force = np.atleast_1d(force)

    # 按比例选取 N 个小于受力阈值的构型索引；若构型索引数小于 N，则返回所有构型索引
    indices = np.where(force < force_threshold)[0]
    if len(indices) < num_selected:
        final_indices = indices.tolist()
    else:
        indices_selected = np.linspace(0, len(indices) - 1, num_selected, dtype=int)
        final_indices = indices[indices_selected].tolist()

    return final_indices


def filter_frames(frames: list, force_threshold: float, num_selected: int = 5):
    """筛选构型"""

    indices = force_exceeds_threshold(force_threshold, num_selected)

    filtered_frames = [frames[i] for i in indices]

    removed_frames = [frames[i] for i in range(len(frames)) if i not in indices]
    removed_frames_indices = [i + 1 for i in range(len(frames)) if i not in indices]

    return filtered_frames, removed_frames, removed_frames_indices


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Select/Filter frames from extxyz file according to the force threshold.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "input_xyz_fn",
        default="train.xyz",
        help="input xyz filename (default: train.xyz)",
    )
    parser.add_argument(
        "output_xyz_fn",
        default="train_filtered.xyz",
        help="output xyz filename (default: train_filtered.xyz)",
    )
    parser.add_argument(
        "removed_xyz_fn",
        default="train_removed.xyz",
        help="removed xyz filename (default: train_removed.xyz)",
    )
    parser.add_argument(
        "-f",
        "--force_threshold",
        type=float,
        default=5.0,
        metavar="FLOAT",
        help="force threshold",
    )
    parser.add_argument(
        "-n",
        "--num_selected",
        type=int,
        default=5,
        metavar="N",
        help="number of frames to select",
    )

    args = parser.parse_args()

    input_xyz_fn = args.input_xyz_fn
    output_xyz_fn = args.output_xyz_fn
    removed_xyz_fn = args.removed_xyz_fn
    force_threshold = args.force_threshold
    num_selected = args.num_selected

    frames = parse_xyz_file(input_xyz_fn)
    filtered_frames, removed_frames, removed_frames_indices = filter_frames(
        frames, force_threshold, num_selected
    )
    write_xyz_file(filtered_frames, output_xyz_fn)
    write_xyz_file(removed_frames, removed_xyz_fn)

    print(f"\nRemoved frame indices (starting from 1):\n")
    print(f"{removed_frames_indices}.")

    print(
        f"\nFiltered structures saved to {output_xyz_fn}, removed structures saved to {removed_xyz_fn}."
    )
