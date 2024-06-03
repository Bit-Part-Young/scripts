"""将 VASP POSCAR Direct 坐标范围限制在 0-1 之间"""

import sys

from ase.io import read, write


def uniform_direct_coords(filename: str):
    """将 VASP POSCAR Direct 坐标范围限制在 0-1 之间"""

    atoms = read(filename)
    scale_positions = atoms.get_scaled_positions()
    atoms.set_scaled_positions(scale_positions)

    write(
        filename,
        images=atoms,
        format="vasp",
        direct=True,
        sort=False,
    )

    print("The direct coordinates range have been uniformed to 0-1.")


if __name__ == "__main__":
    filename = sys.argv[1]

    uniform_direct_coords(filename=filename)
