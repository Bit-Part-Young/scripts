"""
NEB 插值（使用 ASE NEB 模块）

reference: https://github.com/DYang90/ase_extend/blob/main/makeneb.py
"""

import argparse
import os

from ase.atoms import Atoms
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import BFGS, FIRE, LBFGS, MDMin

# [ ] 待完善


def neb_interpolate(
    initial_structure_fn: str = "POSCAR",
    final_structure_fn: str = "POSCAR",
    method: str = "idpp",
    nimage: int = 6,
    spring: float = 0.1,
    fmax: float = 0.1,
    steps: int = 100,
    optimizer: str = "MDMin",
):

    initial_atoms = read(initial_structure_fn)
    final_atoms = read(final_structure_fn)

    images = [initial_atoms] + [initial_atoms.copy() for i in range(nimage - 2)] + [final_atoms]

    neb = NEB(images, k=spring)
    # 可以按照https://wiki.fysik.dtu.dk/ase/ase/neb.html，最终NEB中的image将被更新。
    if method == "linear":
        neb.interpolate(method="linear")
    elif method == "idpp":
        neb.interpolate(method="idpp")

        optimizer_dict = {"MDMin": MDMin, "BFGS": BFGS, "LBFGS": LBFGS, "FIRE": FIRE}
        neb.idpp_interpolate(fmax=fmax, optimizer=optimizer_dict[optimizer], steps=steps)

    return images


def write_guess(images: list[Atoms], output: bool = False):

    for index, image in enumerate(images):
        os.makedirs(f"{index:02d}", exist_ok=True)
        write(f"{index:02d}/POSCAR", image)

    if output:
        write("XDATCAR", images, format="vasp-xdatcar", append=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="ASE NEB Interpolation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--images",
        nargs="+",
        type=str,
        help="Images defining path from initial to final state.",
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=["linear", "idpp"],
        default="linear",
        help="Interpolate method to initial guess. Default: linear",
    )
    parser.add_argument(
        "-n",
        "--nimage",
        type=int,
        default=6,
        help="Number of images in a band. Default: 6",
    )
    parser.add_argument("-o", "--output", action="store_true", help="Whether or not write XDATCAR.")

    parser.add_argument(
        "--nstep",
        type=int,
        default=100,
        help="Max number of iteration steps. Default: 100",
    )
    parser.add_argument(
        "--spring",
        type=float,
        default=0.1,
        help="Fraction of spring force. Default: 0.1",
    )
    parser.add_argument(
        "--fmax",
        type=float,
        default=0.1,
        help="Maximum force of converge thershould. Default: 0.1",
    )
    parser.add_argument(
        "--optimizer",
        type=str,
        choices=["MDMin", "BFGS", "LBFGS", "FIRE"],
        default="MDMin",
        help="Optimizer using in IDPP iteration. Default: MDMin",
    )

    args = parser.parse_args()

    images = neb_interpolate(args)
    write_guess(args, images)
