"""
生成 GPUMD 声子计算的 输入文件：model.xyz、basis.in 和 kpoints.in 文件

reference: https://github.com/tang070205/tools/blob/main/create_phonon_compare.py
"""

import argparse

from ase.io import read, write


def basisin_generation(structure_fn: str = "POSCAR", cellsize: int = 4):
    """生成 GPUMD 声子计算的 basis.in 文件"""

    primitive_structure = read(structure_fn)

    supercell = primitive_structure * (cellsize, cellsize, cellsize)
    write("model.xyz", supercell)

    with open("basis.in", "w") as f:
        f.write(f"{len(primitive_structure)}\n")
        for i, mass in enumerate(primitive_structure.get_masses()):
            f.write(f"{i} {mass}\n")
        for _ in range(cellsize * cellsize * cellsize):
            for i in range(len(primitive_structure)):
                f.write(f"{i}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate GPUMD input files for phonon calculation (model.xyz, basis.in, kpoints.in).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        default="POSCAR",
        help="primitive structure filename",
    )
    parser.add_argument(
        "--cellsize",
        type=int,
        default=4,
        const=4,
        nargs="?",
        metavar="N",
        help="cellsize",
    )

    args = parser.parse_args()
    basisin_generation(args.structure_fn, args.cellsize)
