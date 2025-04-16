#!/usr/bin/env python3

"""NEP 描述符（未归一化） PCA 2 维绘制"""

# [ ] 是否可进一步优化以节约耗时

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from calorine.nep import get_descriptors
from sklearn.decomposition import PCA

np.set_printoptions(suppress=True, precision=10)


def nep_pca_plot(
    descriptors_fn: str = "descriptors.npy",
    output_fn: str = "pca.png",
):
    """NEP 描述符（未归一化） PCA 2 维绘制"""

    descriptors = np.load(descriptors_fn)

    fig, ax = plt.subplots()

    pca = PCA(n_components=2)
    pc = pca.fit_transform(descriptors)

    ax.scatter(
        pc[:, 0],
        pc[:, 1],
        alpha=0.5,
        s=10,
        label="structures",
    )

    ax.set_xlabel(f"PAC dim0 - Var={pca.explained_variance_ratio_[0]:.2f}")
    ax.set_ylabel(f"PAC dim1 - Var={pca.explained_variance_ratio_[1]:.2f}")

    fig.savefig(output_fn)

    print("\nNEP descriptors PCA 2D figure generated!")


def get_nep_descriptors(
    xyz_fn: str = "train.xyz",
    model_fn: str = "nep.txt",
    descriptors_fn: str = "descriptors.npy",
):
    """获取 NEP 描述符（未归一化）"""

    atoms_list = read(xyz_fn, index=":", format="extxyz")

    descriptors_list = []
    for atoms in atoms_list:
        descriptors_list.append(
            get_descriptors(structure=atoms, model_filename=model_fn)
        )

    # 行合并所有描述符
    descriptors_array = np.concatenate(descriptors_list, axis=0)

    print(f"Total number of atoms in dataset: {descriptors_array.shape[0]}")
    print(f"Number of descriptor components:  {descriptors_array.shape[1]}")

    # 数据集很大时，建议保存成文件
    np.save(descriptors_fn, descriptors_array)


def main():
    descriptors_fn = "descriptors.npy"
    if not os.path.exists(descriptors_fn):
        get_nep_descriptors(descriptors_fn=descriptors_fn)

    nep_pca_plot(descriptors_fn=descriptors_fn)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("NEP descriptors PCA 2D plot.")

    main()
