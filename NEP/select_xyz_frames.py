#!/usr/bin/env python3

"""
根据原子受力阈值筛选 extxyz 文件中的构型

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


def force_exceeds_threshold(atom_info: str, force_threshold: float):
    """判断原子受力是否超过阈值"""

    atom_info_list = atom_info.split()
    if len(atom_info_list) == 7:
        force_list = atom_info_list[4:7]
    elif len(atom_info_list) == 8:
        force_list = atom_info_list[5:8]

    # 力矢量
    force_vector = np.array(list(map(float, force_list)))
    # 合力
    force_norm = np.linalg.norm(force_vector)

    return force_norm > force_threshold


def filter_frames(frames: list, force_threshold: float) -> tuple[list, list, list]:
    """筛选构型"""

    filtered_frames = []
    removed_frames = []
    removed_frames_indices = []

    for i in range(len(frames)):
        current_frame = frames[i]
        natoms, frame_info, energy, lattice, atoms_info = current_frame

        # 如果力超过阈值，则删除该帧对应的构型
        if any(
            force_exceeds_threshold(atom_info, force_threshold)
            for atom_info in atoms_info
        ):
            # 记录被删除的构型索引，从 1 开始
            removed_frames_indices.append(i + 1)
            removed_frames.append(current_frame)

            continue

        filtered_frames.append(current_frame)

    return filtered_frames, removed_frames, removed_frames_indices


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Select/Filter frames from extxyz file according to the force threshold.",
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
        "-f",
        "--force_threshold",
        type=float,
        default=50.0,
        metavar="FLOAT",
        help="force threshold",
    )
    args = parser.parse_args()

    input_xyz_fn = args.input_xyz_fn
    output_xyz_fn = args.output_xyz_fn
    force_threshold = args.force_threshold

    frames = parse_xyz(input_xyz_fn)
    filtered_frames, removed_frames, removed_frames_indices = filter_frames(
        frames, force_threshold
    )
    write_xyz(filtered_frames, output_xyz_fn)
    write_xyz(removed_frames, "removed.xyz")

    print(f"\nRemoved frame indices (starting from 1):\n")
    print(f"{removed_frames_indices}.")

    print(
        f"\nFiltered structures saved to {output_xyz_fn}, removed structures saved to removed.xyz."
    )
