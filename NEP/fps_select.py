#!/usr/bin/env python3

"""
通过最远点采样选择结构

reference: https://github.com/bigd4/PyNEP/blob/master/examples/plot_select_structure.py
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read, write
from pynep.calculate import NEP
from pynep.select import FarthestPointSample
from sklearn.decomposition import PCA


# [ ] 是否添加训练集的描述符
def fps_select(
    input_xyz_fn: str,
    model_fn: str = "nep.txt",
    output_xyz_fn: str = "selected.xyz",
    kept_xyz_fn: str = "kept.xyz",
    min_distance: float = 0.05,
):
    """通过 最远点采样 选择结构"""

    if os.path.exists(model_fn):
        calculator = NEP(model_fn)
        print("\n")
        print(calculator)
    else:
        print(f"{model_fn} does not exist! Please check.")
        exit()

    atoms_list = read(input_xyz_fn, index=":", format="extxyz")

    descriptors = np.array(
        [np.mean(calculator.get_descriptor(atoms), axis=0) for atoms in atoms_list]
    )

    sampler = FarthestPointSample(min_distance=min_distance)

    selected_indices = sampler.select(descriptors, [])
    selected_indices.sort()
    selected_indices = [int(i) for i in selected_indices]

    # 保存 FPS 选中的构型
    selected_atoms = [atoms_list[i] for i in selected_indices]
    if os.path.exists(output_xyz_fn):
        os.remove(output_xyz_fn)
    write(output_xyz_fn, selected_atoms, format="extxyz", append=True)

    # 保存 FPS 未选中的构型
    kept_indices = sorted(set(range(len(atoms_list))) - set(selected_indices))
    kept_atoms = [atoms_list[i] for i in kept_indices]
    if os.path.exists(kept_xyz_fn):
        os.remove(kept_xyz_fn)
    write(kept_xyz_fn, kept_atoms, format="extxyz", append=True)

    print(f"Number of selected structures: {len(selected_indices)}.")
    print(f"Slected structures saved to {output_xyz_fn}.")
    print(f"Kept structures saved to {kept_xyz_fn}.")
    print(f"\nIndex of selected structures: {selected_indices}.")

    return descriptors, selected_indices


def plot_pca(descriptors, selected_indices):
    reducer = PCA(n_components=2)
    reducer.fit(descriptors)
    proj = reducer.transform(descriptors)

    selected_proj = reducer.transform(
        np.array([descriptors[i] for i in selected_indices])
    )

    fig, ax = plt.subplots()

    ax.scatter(proj[:, 0], proj[:, 1], label="All")

    ax.scatter(selected_proj[:, 0], selected_proj[:, 1], label="Selected")

    ax.legend()

    fig.savefig("fps_select.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_xyz_fn", type=str)
    parser.add_argument("model_fn", type=str, default="nep.txt")
    parser.add_argument("output_xyz_fn", type=str, default="selected.xyz")
    parser.add_argument("-md", "--min_distance", type=float, default=0.05)
    args = parser.parse_args()

    descriptors, selected_indices = fps_select(
        input_xyz_fn=args.input_xyz_fn,
        model_fn=args.model_fn,
        output_xyz_fn=args.output_xyz_fn,
        min_distance=args.min_distance,
    )
    plot_pca(descriptors, selected_indices)
