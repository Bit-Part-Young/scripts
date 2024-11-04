"""
获取 pymatgen 推荐的赝势文件
注：无法获取 W_sv 赝势文件；POSCAR 中的元素会重新排序
"""

import os
import sys

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import MPMetalRelaxSet


def get_potcar_pymatgen(
    structure_fn: str = "POSCAR",
):
    """获取 pymatgen 推荐的赝势文件"""

    if os.path.exists("POTCAR"):
        os.remove("POTCAR")

    flag = False
    if structure_fn != "POSCAR":
        flag = True

    structure = Structure.from_file(structure_fn)

    settings = MPMetalRelaxSet(structure)
    settings.write_input(".")

    os.remove("INCAR")
    os.remove("KPOINTS")

    if flag:
        os.remove("POSCAR")

    print("\nPOTCAR recommended by pymatgen generated.")


if __name__ == "__main__":

    structure_fn = sys.argv[1] if len(sys.argv) > 1 else "POSCAR"

    get_potcar_pymatgen(
        structure_fn=structure_fn,
    )
