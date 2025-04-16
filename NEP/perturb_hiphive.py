"""使用 hiphive 生成扰动结构（原子位置 扰动）"""

import argparse
import time

from ase.atoms import Atoms
from ase.io import read, write
from hiphive.structure_generation import (
    generate_mc_rattled_structures,
    generate_rattled_structures,
)


def perturb_hiphive(
    atoms: Atoms,
    num_perturb: int = 100,
    rattle_std: float = 0.1,
    d_min: float = 1.5,
    n_iter: int = 10,
    output_fn: str = "perturb_hiphive.xyz",
):
    # seed 设置为当前时间
    seed = int(time.time())

    atoms_list = generate_mc_rattled_structures(
        atoms,
        n_structures=num_perturb,
        rattle_std=rattle_std,
        d_min=d_min,
        seed=seed,
        n_iter=n_iter,
    )

    write(output_fn, atoms_list, format="extxyz", append=True)


if __name__ == "__main__":
    atoms = read("Nb_bcc_dup.vasp", format="vasp")
    perturb_hiphive(atoms)

    print("Work is done!")
