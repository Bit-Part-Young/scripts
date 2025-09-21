#!/usr/bin/env python3

"""给特定的 group 添加权重（主要是 stacking fault, SIA 等缺陷）"""

import argparse


def parse_xyz_file(filename: str) -> list:
    """解析 extxyz 文件"""

    with open(filename, "r") as file:
        lines = file.readlines()

    frames = []
    i = 0
    while i < len(lines):
        num_atoms = int(lines[i].strip())
        frame_info = lines[i + 1].strip()

        group_str = frame_info.split('group="')[1].split('"')[0]

        atoms = lines[i + 2 : i + 2 + num_atoms]
        frames.append((num_atoms, frame_info, group_str, atoms))

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


def add_weight(frames: list, group_info: str, weight: float):
    """添加权重"""

    frames_with_weight = []
    for frame in frames:
        num_atoms, frame_info, group_str, atoms = frame
        if group_str == group_info:
            frame_info = frame_info.replace(
                'tag="train"', f'tag="train" weight={weight}'
            )
        frames_with_weight.append((num_atoms, frame_info, atoms))

    return frames_with_weight


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add weight to specific group in extxyz file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument(
        "input_filename",
        type=str,
        default="train.xyz",
        help="input xyz filename (default: train.xyz)",
    )
    parser.add_argument(
        "output_filename",
        type=str,
        default="train_filtered.xyz",
        help="output xyz filename (default: train_filtered.xyz)",
    )
    parser.add_argument(
        "--group",
        default="stacking fault",
        metavar="STR",
        help="group name (default: stacking fault)",
    )
    parser.add_argument(
        "--weight", type=float, default=2.0, metavar="FLOAT", help="weight"
    )
    args = parser.parse_args()

    input_filename = args.input_filename
    output_filename = args.output_filename
    group = args.group
    weight = args.weight

    frames = parse_xyz_file(input_filename)
    frames_with_weight = add_weight(frames, group_info=group, weight=weight)
    write_xyz_file(frames_with_weight, output_filename)

    print(
        f"Added weight={weight} to group {group} in {input_filename} and saved to {output_filename}."
    )
