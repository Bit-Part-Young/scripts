"""
BCC、FCC、HCP 结构的自间隙原子构型生成

reference: https://github.com/deepmodeling/APEX/blob/main/apex/core/property/Interstitial.py
"""

import argparse

import numpy as np
from ase.atoms import Atoms
from ase.cell import Cell
from ase.io import write

# BCC 自间隙原子位置
bcc_sia_dict = {
    "tetrahedral": [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.25, 0.5, 0],
    ],
    "octahedral": [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0],
    ],
    "crowdion": [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.25, 0.25, 0.25],
    ],
    "<111>dumbbell": [
        [0.0, 0.0, 0.0],
        [1 / 3, 1 / 3, 1 / 3],
        [2 / 3, 2 / 3, 2 / 3],
    ],
    "<110>dumbbell": [
        [0.0, 0.0, 0.0],
        [1 / 4, 3 / 4, 1 / 2],
        [3 / 4, 1 / 4, 1 / 2],
    ],
    "<100>dumbbell": [
        [0.0, 0.0, 0.0],
        [1 / 2, 1 / 2, 1 / 6],
        [1 / 2, 1 / 2, 5 / 6],
    ],
}

# [ ] 待理解各类型原子位置含义
# FCC 自间隙原子位置
fcc_sia_dict = {
    "tetrahedral": [
        [0.75, 0.25, 0.25],
    ],
    "octahedral": [
        [1, 0, 0.5],
    ],
    "crowdion": [
        [1, 0.25, 0.25],
    ],
    "<111>dumbbell": {
        [1 - 0.3 / np.sqrt(3), 1 - 0.3 / np.sqrt(3), 1 - 0.3 / np.sqrt(3)],
        [0.3 / np.sqrt(3), 0.3 / np.sqrt(3), 0.3 / np.sqrt(3)],
    },
    "<110>dumbbell": [
        [1, 0.5 + (0.3 / np.sqrt(2)), 0.5 + (0.3 / np.sqrt(2))],
        [1, 0.5 - (0.3 / np.sqrt(2)), 0.5 - (0.3 / np.sqrt(2))],
    ],
    "<100>dumbbell": [
        [1, 0.2, 0.5],
        [1, 0.8, 0.5],
    ],
}


# [ ] 待完善
def sia_generation(symbols: str, a: float, sia_type: str, output_fn: str):
    """BCC、FCC、HCP 结构的自间隙原子构型生成"""

    if "100" in sia_type:
        sia_type = "<100>dumbbell"
    elif "110" in sia_type:
        sia_type = "<110>dumbbell"
    elif "111" in sia_type:
        sia_type = "<111>dumbbell"

    atoms = Atoms(
        symbols=[symbols] * 3,
        cell=Cell.fromcellpar([a, a, a, 90, 90, 90]),
        scaled_positions=bcc_sia_dict[sia_type],
        pbc=True,
    )

    write(output_fn, atoms, format="vasp", direct=True, sort=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Self interstitial structure generation for BCC, FCC, HCP,",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
