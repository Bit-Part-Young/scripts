#!/usr/bin/env python3

"""解析 extxyz 文件并获取 energy, virial/stress, forces 信息"""


import argparse
import os

import numpy as np

np.set_printoptions(precision=10, suppress=True)


def parse_xyz(input_xyz_fn: str) -> list:
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


def parse_info(
    frames: list, atoms_indices: list[int], output_xyz_fn: str | None = None
):
    """解析 energy, virial/stress, forces 信息"""

    frames_selected = []
    for i, frame in enumerate(frames):
        if i not in atoms_indices:
            continue

        natoms, frame_info, atoms_info = frame
        frame_info: str

        print()
        print("-" * 100)

        print(f"The {i} frame info:\n")

        energy_str = frame_info.split("energy=")[1].split()[0]
        energy = float(energy_str)
        print(f"energy: {round(energy, 5)} eV.\n")

        virial = None
        virial_tag_list = ["virial", "Virial", "virials", "Virials", "stress", "Stress"]
        for virial_tag in virial_tag_list:
            if virial_tag in frame_info:
                virial_str = frame_info.split(f'{virial_tag}="')[1].split('"')[0]
                virial = np.array(list(map(float, virial_str.split()))).reshape(3, 3)

                print(f"{virial_tag}:")
                print(virial)
                print()

                break

        if virial is None:
            print("Virial/Stress: N/A.\n")

        forces = None
        if len(atoms_info[0].split()) > 6:
            forces = np.array(
                [list(map(float, atom_info.split()[-3:])) for atom_info in atoms_info]
            )

            print(f"forces (only show first 10 atoms):")
            if forces.shape[0] > 10:
                print(forces[:10])
            else:
                print(forces)
        else:
            print("Forces: N/A.\n")

        frames_selected.append((natoms, frame_info, atoms_info))

    if output_xyz_fn is not None:
        if os.path.exists(output_xyz_fn):
            os.remove(output_xyz_fn)

        with open(output_xyz_fn, "w") as f:
            for natoms, frame_info, atoms_info in frames_selected:
                f.write(f"{natoms}\n")
                f.write(f"{frame_info}\n")
                for atom_info in atoms_info:
                    f.write(f"{atom_info.strip()}\n")

        print(f"\nSelected frames saved to {output_xyz_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse energy, virial and forces info from extxyz file.",
        epilog="Author: SLY.",
    )

    parser.add_argument("input_xyz_fn", help="input xyz filename")
    parser.add_argument(
        "-ai", "--atoms_indices", nargs="+", type=int, metavar="N", help="atoms indices"
    )
    parser.add_argument(
        "-o", "--output_xyz_fn", metavar="FILE", help="output xyz filename"
    )

    args = parser.parse_args()

    frames = parse_xyz(args.input_xyz_fn)
    parse_info(frames, args.atoms_indices, args.output_xyz_fn)
