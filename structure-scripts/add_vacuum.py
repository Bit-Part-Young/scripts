"""
添加真空层
包括：z 方向顶部、底部、底部顶部两端
"""

from ase.atoms import Atoms
from ase.build import bulk

# TODO: 添加底部真空层、底部顶部两端真空层函数


def add_top_vacuum(atoms: Atoms, vacuum: float) -> Atoms:
    """在顶部添加真空层"""

    cell = atoms.get_cell()
    cell[2, 2] += vacuum

    atoms.set_cell(cell)

    return atoms


if __name__ == "__main__":

    atoms = bulk("Al", a=4.039, cubic=True)
    print(atoms.get_cell())
    print(atoms.get_scaled_positions())

    new_atoms = add_top_vacuum(atoms, 10)
    print(new_atoms.get_cell())
    print(new_atoms.get_scaled_positions())
