#!/usr/bin/env python3

"""
根据索引列表筛选构型

reference: https://github.com/brucefan1983/GPUMD/blob/master/tools/Analysis_and_Processing/select_xyz_frames/select_xyz_frames.py
"""

import argparse

import numpy as np


def parse_xyz(input_xyz_fn: str) -> list:
    """解析 extxyz 文件"""

    with open(input_xyz_fn, "r") as file:
        lines = file.readlines()

    frames = []
    i = 0
    while i < len(lines):
        natoms = int(lines[i].strip())
        frame_info = lines[i + 1].strip()

        energy_str = frame_info.split("energy=")[1].split()[0]
        energy = float(energy_str)

        lattice_str = frame_info.split('Lattice="')[1].split('"')[0]
        lattice = np.array(list(map(float, lattice_str.split()))).reshape(3, 3)

        atoms_info = lines[i + 2 : i + 2 + natoms]
        frames.append((natoms, frame_info, energy, lattice, atoms_info))
        i += 2 + natoms
    return frames


def write_xyz(frames: list, output_xyz_fn: str):
    """写入 extxyz 文件"""

    with open(output_xyz_fn, "w") as f:
        for natoms, frame_info, energy, lattice, atoms_info in frames:
            f.write(f"{natoms}\n")
            f.write(f"{frame_info}\n")
            for atom_info in atoms_info:
                f.write(f"{atom_info.strip()}\n")


def filter_frames(frames: list, selected_indices: list) -> tuple[list, list, list]:
    """根据给定索引列表筛选构型"""

    filtered_frames = []
    removed_frames = []
    removed_frames_indices = []

    for i in range(len(frames)):
        current_frame = frames[i]

        # 若索引不在 selected_indices 中，则删除该帧对应的构型
        if i not in selected_indices:
            removed_frames_indices.append(i)
            removed_frames.append(current_frame)

            continue

        filtered_frames.append(current_frame)

    return filtered_frames, removed_frames, removed_frames_indices


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Select/Filter configuration frames by selected indices from extxyz file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        "-si",
        "--selected_indices",
        type=int,
        nargs="+",
        metavar="N",
        help="selected indices (starting from 1)",
    )

    args = parser.parse_args()

    input_xyz_fn = args.input_xyz_fn
    output_xyz_fn = args.output_xyz_fn
    selected_indices = args.selected_indices

    frames = parse_xyz(input_xyz_fn)
    filtered_frames, removed_frames, removed_frames_indices = filter_frames(
        frames, selected_indices
    )
    write_xyz(filtered_frames, output_xyz_fn)
    write_xyz(removed_frames, "removed.xyz")

    print(f"\nRemoved frame indices:\n")
    print(f"{removed_frames_indices}.")

    print(
        f"\nFiltered structures saved to {output_xyz_fn}, removed structures saved to removed.xyz."
    )
