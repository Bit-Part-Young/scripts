#!/usr/bin/env python3

"""按比例选取指定个数的小于受力阈值的构型（用于缺陷、层错构型弛豫轨迹构型的筛选）"""

import argparse

import numpy as np


def parse_xyz_file(filename: str) -> list:
    """解析 extxyz 文件"""

    with open(filename, "r") as file:
        lines = file.readlines()

    frames = []
    i = 0
    while i < len(lines):
        num_atoms = int(lines[i].strip())
        frame_info = lines[i + 1].strip()

        atoms = lines[i + 2 : i + 2 + num_atoms]
        frames.append((num_atoms, frame_info, atoms))
        i += 2 + num_atoms

    return frames


def write_xyz_file(frames: list, output_filename: str):
    """写入 extxyz 文件"""

    with open(output_filename, "w") as file:
        for num_atoms, frame_info, atoms in frames:
            file.write(f"{num_atoms}\n")
            file.write(f"{frame_info}\n")
            for atom in atoms:
                file.write(f"{atom.strip()}\n")


def force_exceeds_threshold(threshold: float = 5.0, num_selected: int = 5) -> list:
    """按比例选取指定个数的小于受力阈值的构型索引"""

    force = np.loadtxt("force_info.dat", skiprows=1)

    # 按比例选取 N 个小于受力阈值的构型索引；若构型索引数小于 N，则返回所有构型索引
    index = np.where(force < threshold)[0]
    if len(index) < num_selected:
        final_index = index.tolist()
    else:
        index_selected = np.linspace(0, len(index) - 1, num_selected, dtype=int)
        final_index = index[index_selected].tolist()

    return final_index


def filter_frames(frames: list, force_threshold: float, num_selected: int = 5):
    """筛选构型"""

    index = force_exceeds_threshold(force_threshold, num_selected)

    filtered_frames = [frames[i] for i in index]

    removed_frames = [frames[i] for i in range(len(frames)) if i not in index]
    removed_frames_indices = [i + 1 for i in range(len(frames)) if i not in index]

    return filtered_frames, removed_frames, removed_frames_indices


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Select/Filter frames from exyz configuration file according to the force threshold.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "input_filename",
        type=str,
        default="train.xyz",
        help="input filename (default: train.xyz)",
    )
    parser.add_argument(
        "output_filename",
        type=str,
        default="train_filtered.xyz",
        help="output filename (default: train_filtered.xyz)",
    )
    parser.add_argument(
        "removed_filename",
        type=str,
        default="train_removed.xyz",
        help="removed filename (default: train_removed.xyz)",
    )
    parser.add_argument(
        "-f",
        "--force_threshold",
        type=float,
        default=5.0,
        required=True,
        help="force threshold",
    )
    parser.add_argument(
        "-n",
        "--num_selected",
        type=int,
        default=5,
        help="number of frames to select",
    )

    args = parser.parse_args()

    input_filename = args.input_filename
    output_filename = args.output_filename
    removed_filename = args.removed_filename
    force_threshold = args.force_threshold
    num_selected = args.num_selected

    frames = parse_xyz_file(input_filename)
    filtered_frames, removed_frames, removed_frames_indices = filter_frames(
        frames, force_threshold, num_selected
    )
    write_xyz_file(filtered_frames, output_filename)
    write_xyz_file(removed_frames, removed_filename)

    print(f"\nRemoved frame indices (starting from 1): {removed_frames_indices}.")

    print(
        f"\nFiltered structures saved to {output_filename}, removed structures saved to {removed_filename}."
    )
