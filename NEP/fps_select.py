"""
通过最远点采样选择结构

reference: https://github.com/bigd4/PyNEP/blob/master/examples/plot_select_structure.py
"""

import argparse

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
    min_distance: float = 0.05,
):
    """通过 最远点采样 选择结构"""

    atoms_list = read(input_xyz_fn, index=":", format="extxyz")

    calculator = NEP(model_fn)
    print(calculator)

    descriptors = np.array(
        [np.mean(calculator.get_descriptor(atoms), axis=0) for atoms in atoms_list]
    )

    sampler = FarthestPointSample(min_distance=min_distance)

    selected_indices = sampler.select(descriptors, [])
    selected_atoms = [atoms_list[i] for i in selected_indices]

    write(output_xyz_fn, selected_atoms, format="extxyz", append=True)

    print(f"Number of selected structures: {len(selected_indices)}.")

    return descriptors, selected_indices


def plot_pca(descriptors, selected_indices):
    reducer = PCA(n_components=2)
    reducer.fit(descriptors)
    proj = reducer.transform(descriptors)

    selected_proj = reducer.transform(
        np.array([descriptors[i] for i in selected_indices])
    )

    fig, ax = plt.subplots()

    ax.scatter(proj[:, 0], proj[:, 1], label="all data")

    ax.scatter(selected_proj[:, 0], selected_proj[:, 1], label="selected data")

    ax.legend()

    fig.savefig("select.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_xyz_fn", type=str)
    parser.add_argument("model_fn", type=str, default="nep.txt")
    parser.add_argument("output_xyz_fn", type=str, default="selected.xyz")
    parser.add_argument("-md", "--min_distance", type=float, default=0.05)
    args = parser.parse_args()

    descriptors, selected_indices = fps_select(
        args.input_xyz_fn, args.model_fn, args.output_xyz_fn
    )
    plot_pca(descriptors, selected_indices)
