"""
对原子位置进行扰动、 施加应变使结构发生变形

reference: https://github.com/aboys-cb/NepTrain/blob/master/src/NepTrain/core/perturb/run.py
"""

import numpy as np
from ase.atoms import Atoms
from ase.build import bulk


def perturb_position(
    prim: Atoms,
    min_distance: float,
):
    """对原子位置进行扰动"""

    atoms = prim.copy()
    positions = atoms.get_positions()

    # 添加原子位置的随机微扰
    perturbed_positions = positions + np.random.uniform(
        low=-min_distance,
        high=min_distance,
        size=positions.shape,
    )

    # 更新原子位置
    atoms.set_positions(perturbed_positions)

    return atoms


def generate_strained_structure(
    prim: Atoms,
    strain_lim: list[float],
    min_distance: float,
):
    """正应变"""

    strains = np.random.uniform(*strain_lim, (3,))
    atoms = prim.copy()
    cell_new = prim.cell[:] * (1 + strains)
    atoms.set_cell(cell_new, scale_atoms=True)

    return perturb_position(atoms, min_distance)


def generate_deformed_structure(
    prim: Atoms,
    strain_lim: list[float],
    min_distance: float,
):
    """正应变 + 切应变"""

    R = np.random.uniform(*strain_lim, (3, 3))
    M = np.eye(3) + R

    atoms = prim.copy()
    cell_new = M @ atoms.cell[:]
    atoms.set_cell(cell_new, scale_atoms=True)

    return perturb_position(atoms, min_distance)


if __name__ == "__main__":
    prim = bulk("Al", cubic=True)
    strain_lim = (-0.05, 0.05)

    # 不用科学计数法表示
    np.set_printoptions(suppress=True)

    generate_strained_structure(prim, strain_lim, 0.1)

    generate_deformed_structure(prim, strain_lim, 0.1)
