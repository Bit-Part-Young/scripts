import matplotlib.pyplot as plt
import numpy as np
from ase.io import read, write
from calorine.nep import get_descriptors
from sklearn.decomposition import PCA

descriptors = []

xyz_fn = ""
atoms_list = read(xyz_fn, index=":", format="extxyz")

# Generate a few different rattled structures
for atoms in atoms_list[:50]:
    descriptors.append(get_descriptors(atoms))

# Concatenate all descriptors
all_descriptors = np.concatenate(descriptors, axis=0)
print(f"Total number of atoms in dataset: {all_descriptors.shape[0]}")
print(f"Number of descriptor components:  {all_descriptors.shape[1]}")

# Perform PCA
pca = PCA(n_components=2)
pc = pca.fit_transform(all_descriptors)

p0 = pca.explained_variance_ratio_[0]
p1 = pca.explained_variance_ratio_[1]
print(f"Explained variance for component 0: {p0:.2f}")
print(f"Explained variance for component 1: {p1:.2f}")


# PCA 绘制
# fig, ax = plt.subplots(figsize=(4, 3), dpi=140)

# ax.scatter(pc[:100, 0], pc[:100, 1], alpha=0.5, label="rattled structures")
# ax.scatter(pc[100:, 0], pc[100:, 1], alpha=0.5, label="rattled EV curves")
# ax.set_xlabel("PCA dimension 0")
# ax.set_ylabel("PCA dimension 1")
# ax.legend(frameon=False)

# plt.tight_layout()


# 描述符归一化
descriptor_mean = all_descriptors.mean(axis=0)
descriptor_std = all_descriptors.std(axis=0)
normalized_descriptors = [(d - descriptor_mean) / descriptor_std for d in descriptors]

fig, axes = plt.subplots(ncols=2, figsize=(7, 3), dpi=140)

for k, ax in enumerate(axes):
    pca = PCA(n_components=2)
    if k == 0:
        pc = pca.fit_transform(all_descriptors)
        title = "unnormalized"
    else:
        pc = pca.fit_transform(np.concatenate(normalized_descriptors, axis=0))
        title = "normalized"
    ax.scatter(pc[:100, 0], pc[:100, 1], alpha=0.5, label="rattled structures")
    ax.scatter(pc[100:, 0], pc[100:, 1], alpha=0.5, label="rattled EV curves")
    ax.set_xlabel(f"PCA dimension 0 - Var={pca.explained_variance_ratio_[0]:.2f}")
    ax.set_ylabel(f"PCA dimension 1 - Var={pca.explained_variance_ratio_[1]:.2f}")
    ax.set_title(title)
ax.legend(frameon=False)

plt.tight_layout()
fig.align_ylabels()

