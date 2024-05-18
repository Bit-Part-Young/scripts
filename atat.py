"""实现 ATAT str.out 文件读取"""

import io

import numpy as np
from ase import Atoms
from ase.io.vasp import read_vasp, write_vasp
from ase.utils import reader


@reader
def read_atat_out(file: str, tot_natoms_threshold: int = 1000):
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


# TODO 实现枚举的 str.out 文件批量读取
@reader
def read_atat_enum_out(file: str, index=-1):
    fd = file

    images = fd.read().split("end\n").remove("\n")

    if isinstance(index, int):
        atoms = read_atat_out(io.StringIO(images[index]))

    return atoms
