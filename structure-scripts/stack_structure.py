"""
在 z 方向上堆叠/合并两个相同的结构

reference: https://github.com/wangchr1617/learning/blob/main/scripts/model/combine_vasp_structures.py
"""

import numpy as np
from ase.io import read, write


def stack_structure(
    structure1_fn: str,
    structure2_fn: str,
    spacing: float = 2.0,
    output_fn: str = "AB.vasp",
):
    """在 z 方向上堆叠/合并两个相同的结构（第一个结构在下层，第二个结构在上层）"""

    atoms1 = read(structure1_fn)
    atoms2 = read(structure2_fn)

    cell1 = atoms1.cell
    cell2 = atoms2.cell
    z_cell1 = cell1[2, 2]

    cell_stack = np.zeros((3, 3))
    cell_stack[0:2] = cell1[0:2]
    cell_stack[2] = cell1[2] + cell2[2] + [0, 0, spacing]

    atoms1.set_cell(cell_stack, scale_atoms=False)
    atoms2.set_cell(cell_stack, scale_atoms=False)

    # TODO: 在第一个结构的 z 方向晶格常数 还是 z 方向的最大坐标值
    z_translation = z_cell1 - atoms2.positions[:, 2].min() + spacing
    # z_translation = atoms1.positions[:, 2].max() - atoms2.positions[:, 2].min() + spacing

    atoms2.translate((0, 0, z_translation))

    atoms_stack = atoms1 + atoms2

    write(
        output_fn,
        atoms_stack,
        format="vasp",
        direct=True,
        sort=True,
    )


if __name__ == "__main__":

    structure1_fn = "POSCAR"
    structure2_fn = "POSCAR"

    stack_structure(structure1_fn, structure2_fn)
