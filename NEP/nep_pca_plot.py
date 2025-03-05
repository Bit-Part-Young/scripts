"""NEP 未归一化 & 归一化 描述符 PCA 2 维绘制

测试 8762 frames, 509917 atoms 的构型训练集 xyz 文件
Master 运行耗时 ~4 min（其中 描述符获取、处理耗时约 ~2 min，PCA 绘图耗时约 ~2 min）
"""

# [ ] 是否可进一步优化以节约耗时

import os

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from calorine.nep import get_descriptors
from matplotlib.axes import Axes
from sklearn.decomposition import PCA
from spt.plot_params import set_roman_plot_params

# np.set_printoptions(suppress=True, precision=10)


def nep_pca_plot(
    all_descriptors_fn: str = "all_descriptors.npy",
    all_normalized_descriptors_fn: str = "all_normalized_descriptors.npy",
):
    """NEP 未归一化 & 归一化 描述符 PCA 2 维绘制"""

    all_descriptors = np.load(all_descriptors_fn)
    all_normalized_descriptors = np.load(all_normalized_descriptors_fn)

    set_roman_plot_params()
    fig, axes = plt.subplots(ncols=2, figsize=(12, 6))

    for i, ax in enumerate(axes):
        ax: Axes
        pca = PCA(n_components=2)
        if i == 0:
            pc = pca.fit_transform(all_descriptors)
            title = "Unnormalized"
        else:
            pc = pca.fit_transform(all_normalized_descriptors)
            title = "Normalized"

        ax.scatter(pc[:, 0], pc[:, 1], alpha=0.5, label="structures")

        ax.set_title(title)
        ax.set_xlabel(f"PAC dim0 - Var={pca.explained_variance_ratio_[0]:.2f}")
        ax.set_ylabel(f"PAC dim1 - Var={pca.explained_variance_ratio_[1]:.2f}")

        # ax.legend(frameon=False)

    plt.tight_layout()
    fig.align_ylabels()

    fig.savefig("pca_plot.png")

    print("PCA 2D unnormalized & normalized figure generated!")


def get_nep_descriptors(
    xyz_fn: str = "train.xyz",
    model_fn: str = "nep.txt",
):
    """获取 NEP 未归一化 & 归一化 描述符"""

    atoms_list = read(xyz_fn, index=":", format="extxyz")

    descriptors = []
    for atoms in atoms_list:
        descriptors.append(get_descriptors(structure=atoms, model_filename=model_fn))

    # 行合并所有描述符
    all_descriptors = np.concatenate(descriptors, axis=0)

    print(f"Total number of atoms in dataset: {all_descriptors.shape[0]}")
    print(f"Number of descriptor components:  {all_descriptors.shape[1]}")

    # 描述符归一化
    descriptor_mean = all_descriptors.mean(axis=0)
    descriptor_std = all_descriptors.std(axis=0)
    normalized_descriptors = [(d - descriptor_mean) / descriptor_std for d in descriptors]
    all_normalized_descriptors = np.concatenate(normalized_descriptors, axis=0)

    # 数据集很大时，建议保存成文件
    np.save("all_descriptors.npy", all_descriptors)
    np.save("all_normalized_descriptors.npy", all_normalized_descriptors)


def main():

    if not os.path.exists("all_descriptors.npy"):
        get_nep_descriptors()

    nep_pca_plot()


if __name__ == "__main__":
    main()
