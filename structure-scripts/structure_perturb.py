"""
生成具有随机晶格扰动（晶格长度 + 晶格角度）的结构

reference: https://github.com/wangchr1617/learning/blob/main/scripts/model/generate_random_perturbed_structures.py
"""

import os

import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar


def generate_random_perturbed_structures(
    structure_fn: str = "POSCAR",
    tolerance=None | list[int],
    nstructure: int = 5,
):
    """生成具有随机晶格扰动（晶格长度 + 晶格角度）的结构"""

    if tolerance is None:
        tolerance = [3, 3]

    origin_structure = Poscar.from_file(structure_fn, check_for_POTCAR=False).structure
    # origin_structure = Structure.from_prototype(prototype="fcc", species=["Al"], a=4.05)

    deformed_structure = origin_structure.copy()

    lengths = origin_structure.lattice.abc
    angles = origin_structure.lattice.angles
    origin_lattice = np.array(lengths + angles)

    delta, sigma = map(int, tolerance)
    Dlength = range(-delta, delta + 1)
    Dangle = range(-sigma, sigma + 1)

    Pfile = "perturb"
    os.makedirs(Pfile, exist_ok=True)

    # 生成具有随机扰动的结构
    for rnd in range(nstructure):
        # 随机选择长度扰动
        lengths_loss = np.random.choice(Dlength, 3) / 100
        lengths_loss = np.array(lengths) * lengths_loss
        # 随机选择角度扰动
        angle_loss = np.random.choice(Dangle, 3)
        loss = np.append(lengths_loss, angle_loss)
        Perturbed_lattice_array = np.array(origin_lattice) - loss

        deformed_structure.lattice = Lattice.from_parameters(*Perturbed_lattice_array)

        out_poscar = Poscar(
            deformed_structure,
            comment=f"{rnd+1}_Perturbed structure",
        )
        out_poscar.write_file(f"{Pfile}/perturbed_{rnd+1}.vasp")

    print(f"{nstructure} structures with random perturbation in lattice were generated! Bye!")


if __name__ == "__main__":
    structure_fn = "POSCAR"
    tolerance = [3, 3]
    nstructure = 5

    generate_random_perturbed_structures(
        structure_fn=structure_fn,
        tolerance=tolerance,
        nstructure=nstructure,
    )
