#!/usr/bin/env python3

"""获取构型（单帧或多帧轨迹）中原子对的最小和最大距离"""

import argparse

import numpy as np
from ase.atoms import Atoms
from ase.io import read


def read_frames(structure_fn: str) -> list[Atoms]:
    """读取全部帧（LAMMPS data 需显式指定格式，其余由 ASE 自动识别）"""

    if structure_fn.rsplit(".", 1)[-1] in ("lmp", "data", "lammps-data"):
        return [read(structure_fn, format="lammps-data", sort_by_id=False)]
    return read(structure_fn, index=":")


def get_distance(atoms: Atoms, mic: bool = False) -> tuple[float, float]:
    """返回原子对 (i < j) 的最小和最大距离（重叠原子的距离 0 也会被计入）"""

    distances = atoms.get_all_distances(mic=mic)[np.triu_indices(len(atoms), k=1)]

    return float(distances.min()), float(distances.max())


def main(structure_fn: str = "POSCAR", mic: bool = False):
    frames = read_frames(structure_fn)
    valid = [k for k, atoms in enumerate(frames) if len(atoms) > 1]

    if not valid:
        print("No frame with more than one atom, no atomic pair.")
        return

    if len(frames) == 1:
        min_dist, max_dist = get_distance(frames[0], mic=mic)
        print(f"Minimum distance: {min_dist:.3f} Å.")
        print(f"Maximum distance: {max_dist:.3f} Å.")
        return

    dists = np.array([get_distance(frames[k], mic=mic) for k in valid]).round(3)
    print(f"\nNumber of frames: {len(frames)}; with more than one atom: {len(valid)}.\n")
    print("Minimum distance list (Å):")
    print(dists[:, 0].tolist())
    print("Maximum distance list (Å):")
    print(dists[:, 1].tolist())
    print(
        f"Global minimum distance: {dists[:, 0].min():.3f} Å "
        f"(frame {valid[dists[:, 0].argmin()]}). Average: {dists[:, 0].mean():.3f} Å."
    )
    print(
        f"Global maximum distance: {dists[:, 1].max():.3f} Å "
        f"(frame {valid[dists[:, 1].argmax()]}, 0-based)."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get the minimum and maximum atomic pair distances of a structure "
        "or of every frame in a trajectory (extxyz, XDATCAR, OUTCAR, ...).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument("structure_fn", nargs="?", default="POSCAR", help="structure filename")
    parser.add_argument(
        "-mic", "--mic", action="store_true", help="use minimum image convention (periodic images)"
    )
    args = parser.parse_args()

    main(args.structure_fn, mic=args.mic)
