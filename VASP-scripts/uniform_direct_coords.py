"""将 VASP POSCAR 中原子坐标分数范围限制在 0-1 之间"""

import sys

from ase.io.vasp import read_vasp, write_vasp


def uniform_direct_coords(filename: str):
    """将 VASP POSCAR 中原子坐标分数范围限制在 0-1 之间"""

    atoms = read_vasp(filename)
    # scale_positions 坐标值范围在 0-1 之间
    scale_positions = atoms.get_scaled_positions()

    atoms.set_scaled_positions(scale_positions)

    write_vasp(filename, atoms, direct=True, sort=False)

    print("The direct coordinates range have been uniformed to 0-1.")


if __name__ == "__main__":
    filename = sys.argv[1]

    uniform_direct_coords(filename=filename)
