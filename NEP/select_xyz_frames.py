#!/usr/bin/env python3

"""
根据力阈值筛选 extxyz 文件中的构型

reference: https://github.com/brucefan1983/GPUMD/blob/master/tools/Analysis_and_Processing/select_xyz_frames/select_xyz_frames.py
"""

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

        energy_str = frame_info.split("energy=")[1].split()[0]
        energy = float(energy_str)

        lattice_str = frame_info.split('Lattice="')[1].split('"')[0]
        lattice = np.array(list(map(float, lattice_str.split()))).reshape(3, 3)

        atoms = lines[i + 2 : i + 2 + num_atoms]
        frames.append((num_atoms, frame_info, energy, lattice, atoms))
        i += 2 + num_atoms
    return frames


def write_xyz_file(frames: list, output_filename: str):
    """写入 extxyz 文件"""

    with open(output_filename, "w") as file:
        for num_atoms, frame_info, energy, lattice, atoms in frames:
            file.write(f"{num_atoms}\n")
            file.write(f"{frame_info}\n")
            for atom in atoms:
                file.write(f"{atom.strip()}\n")


def force_exceeds_threshold(atom_line: str, threshold: float) -> bool:
    """判断原子受力是否超过阈值"""

    if threshold == "not":
        return False

    forces = np.array(list(map(float, atom_line.split()[4:7])))
    # 合力
    total_force = np.linalg.norm(forces)

    return total_force > threshold


def filter_frames(frames: list, force_threshold: float) -> tuple[list, list]:
    """筛选构型"""

    filtered_frames = []
    removed_frames = []
    removed_frames_indices = []

    for i in range(len(frames)):
        current_frame = frames[i]
        num_atoms, frame_info, energy, lattice, atoms = current_frame

        # 如果力超过阈值，则删除该帧对应的构型
        if force_threshold != "not" and any(
            force_exceeds_threshold(atom, force_threshold) for atom in atoms
        ):
            # 记录被删除的构型索引，从 1 开始
            removed_frames_indices.append(i + 1)
            removed_frames.append(current_frame)
            continue

        filtered_frames.append(current_frame)

    return filtered_frames, removed_frames, removed_frames_indices


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Select/Filter frames from exyz configuration file according to the force threshold.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        "-f",
        "--force_threshold",
        type=float,
        default=50.0,
        required=True,
        help="force threshold",
    )
    args = parser.parse_args()

    input_filename = args.input_filename
    output_filename = args.output_filename
    force_threshold = args.force_threshold

    frames = parse_xyz_file(input_filename)
    filtered_frames, removed_frames, removed_frames_indices = filter_frames(
        frames, force_threshold
    )
    write_xyz_file(filtered_frames, output_filename)
    write_xyz_file(removed_frames, "removed.xyz")

    print(f"\nRemoved frame indices (starting from 1): {removed_frames_indices}")

    print(
        f"\nFiltered structures saved to {output_filename}, removed structures saved to removed.xyz"
    )
