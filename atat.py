"""ATAT str.out 文件的构型读取并转换为 ase Atoms 对象"""

import os
import re
from typing import List, Union

import numpy as np
from ase import Atoms
from ase.io import string2index
from ase.utils import reader


@reader
def read_atat_out(
    file: str,
    tot_natoms_threshold: int = 1000,
) -> Atoms:
    """读取 str.out 文件中的构型并转换为 ase Atoms 对象"""

    fd = file

    coordinate_vectors = []
    for _ in range(3):
        s = fd.readline().split()
        floatvect = float(s[0]), float(s[1]), float(s[2])
        coordinate_vectors.append(floatvect)

    lattice_vectors = []
    for _ in range(3):
        s = fd.readline().split()
        floatvect = float(s[0]), float(s[1]), float(s[2])
        lattice_vectors.append(floatvect)

    basis_vectors = np.array(coordinate_vectors) @ np.array(lattice_vectors)

    atoms_pos = []
    atom_symbols = []
    for _ in range(tot_natoms_threshold):
        ac = fd.readline().split()
        if not ac:
            break
        atoms_pos.append([float(ac[0]), float(ac[1]), float(ac[2])])
        atom_symbols.append(ac[3])
    atoms_pos = np.array(atoms_pos)

    cartesian = False
    if np.max(np.abs(atoms_pos)) > 1.0:
        cartesian = True

    atoms = Atoms(symbols=atom_symbols, cell=basis_vectors, pbc=True)
    if cartesian:
        atoms.set_positions(atoms_pos)
    else:
        atoms.set_scaled_positions(atoms_pos)

    return atoms


def read_atat_enum_out(
    file: str,
    index: Union[int, slice, str] = -1,
) -> Union[Atoms, List[Atoms]]:
    """读取枚举的 str.out 文件构型，并实现 `ase.io.read()` 的构型索引功能"""

    with open(file, "r") as fs:
        text = fs.read()

    # 匹配每个 end 前的内容，忽略空白行
    pattern = r"(?:(?!\n\n).)*?\nend"
    matches = re.findall(pattern, text, re.DOTALL)

    if isinstance(index, str):
        try:
            index = string2index(index)
        except ValueError:
            pass

    def pattern2atoms(match: str) -> Atoms:
        match = re.sub(r"\nend", "", match).strip()

        # 将匹配的内容写入临时文件，最后删除
        tmp_fn = "tmp.out"
        with open(tmp_fn, "w") as f:
            f.write(match)

        atoms = read_atat_out(tmp_fn)

        os.remove(tmp_fn)

        return atoms

    if isinstance(index, (slice, str)):
        atoms_list = []
        for match in matches[index]:
            atoms = pattern2atoms(match)
            atoms_list.append(atoms)

        return atoms_list
    else:
        match = matches[index]
        atoms = pattern2atoms(match)

        return atoms


if __name__ == "__main__":
    fn = "str_enum.out"

    atoms_list = read_atat_enum_out(fn, index="4:6")
    print(atoms_list)
