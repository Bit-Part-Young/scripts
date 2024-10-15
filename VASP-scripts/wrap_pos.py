"""将 VASP POSCAR 原子 Direct 坐标范围 wrap 在 0-1 之间"""

import sys

from ase.io import read, write


def wrap_pos(filename: str):
    """将 VASP POSCAR 原子 Direct 坐标 wrap 在 0-1 之间"""

    atoms = read(filename)
    atoms.wrap()

    write(
        filename,
        images=atoms,
        format="vasp",
        direct=True,
        sort=False,
    )

    print("The position has been wrap to 0-1.")


if __name__ == "__main__":
    filename = sys.argv[1]

    wrap_pos(filename=filename)
