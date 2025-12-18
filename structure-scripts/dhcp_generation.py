#!/usr/bin/env python3

"""Generate DHCP structure."""

import argparse

import numpy as np
from ase.atoms import Atoms
from ase.io import write


def dhcp_generation(symbol: str, a: float, c: float):
    """生成 DHCP 结构"""

    cell = ((a, 0, 0), (-0.5 * a, 0.5 * np.sqrt(3) * a, 0), (0, 0, 2 * c))
    scaled_positions = [
        (0 / 3, 0 / 3, 0.0),
        (1 / 3, 2 / 3, 0.25),
        (0 / 3, 0 / 3, 0.5),
        (2 / 3, 1 / 3, 0.75),
    ]
    atoms = Atoms(4 * symbol, cell=cell, scaled_positions=scaled_positions)

    write(f"POSCAR", atoms, format="vasp", direct=True)

    return atoms


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Generate DHCP structure.")

    parser.add_argument("-e", "--symbol", help="element symbol")
    parser.add_argument(
        "-lc", "--lattice-constant", type=float, nargs=2, help="lattice constant"
    )

    args = parser.parse_args()

    dhcp_generation(args.symbol, args.lattice_constant[0], args.lattice_constant[1])
