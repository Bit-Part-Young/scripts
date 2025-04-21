#!/usr/bin/env python3

"""
生成应变结构

reference: http://link.aps.org/ supplemental/10.1103/PhysRevMaterials.7.103602
"""


import argparse

import numpy as np
from ase.io import write
from pymatgen.analysis.elasticity.strain import Deformation, Strain
from pymatgen.core.structure import Structure


def get_default_strain_states(order=2):
    """生成 'strain-states' list"""

    inds = [(i,) for i in range(6)]
    # Note that these strain states may not be minimal
    if order > 2:
        inds.extend([(0, i) for i in range(1, 5)] + [(1, 2), (3, 4), (3, 5), (4, 5)])
        if order > 3:
            inds.extend(
                [
                    (0, 1, 2),
                    (0, 1, 3),
                    (0, 1, 4),
                    (0, 1, 5),
                    (0, 2, 3),
                    (0, 2, 4),
                    (0, 2, 5),
                    (1, 2, 3),
                    (1, 2, 4),
                    (1, 2, 5),
                    (2, 3, 4),
                    (2, 3, 5),
                    (2, 4, 5),
                    (3, 4, 5),
                ]
            )
            if order > 4:
                raise ValueError(
                    "Standard deformations for tensors higher than rank 4 not yet determined"
                )

    strain_states = np.zeros((len(inds), 6))
    for n, i in enumerate(inds):
        np.put(strain_states[n], i, 1)

    # order=2 时，strain_states 为对角线元素为 1 1 1 2 2 2 的 6x6 矩阵
    # order=3 时，strain_states 14x6 矩阵
    # order=3 时，strain_states 28x6 矩阵
    strain_states[:, 3:] *= 2

    return strain_states.tolist()


def get_deformations(order: int = 2, strain_limit: float = 0.1):
    """生成 deformations list"""

    strains = []
    strain_states = get_default_strain_states(order)

    num_stencils = int((strain_limit * 2) / 0.01) + 1 + (order - 2) * 2
    stencils = [
        np.linspace(
            round(-strain_limit, 2),
            round(strain_limit, 2),
            num_stencils,
        )
    ] * len(strain_states)

    # np.array(stencils).ndim 为 2
    if np.array(stencils).ndim == 1:
        stencils = [stencils] * len(strain_states)

    # strain 共含 6*21=126 个元素；strain 为 3x3 矩阵
    for state, stencil in zip(strain_states, stencils):
        strains.extend([Strain.from_voigt(s * np.array(state)) for s in stencil])

    # 删除 为 0 的 strain
    strains = [strain for strain in strains if not np.isclose(strain, 0).all()]
    # vstrains = [strain.voigt for strain in strains]
    deformations = [s.get_deformation_matrix() for s in strains]

    return deformations


def strained_structure_generation(
    order: int = 2,
    input_fn: str = "POSCAR",
    output_fn: str = "strained.xyz",
    strain_limit: float = 0.1,
):
    """生成应变结构"""

    structure = Structure.from_file(input_fn)

    deformations = get_deformations(order, strain_limit)

    strained_structure_list = [
        defo.apply_to_structure(structure) for defo in deformations
    ]

    for strained_structure in strained_structure_list:
        strained_structure: Structure
        atoms = strained_structure.to_ase_atoms()

        write(output_fn, atoms, format="extxyz", append=True)

    print(f"Total {len(strained_structure_list)} strained structures generated.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate strained structures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "order",
        type=int,
        nargs="?",
        const=2,
        default=2,
        choices=[2, 3, 4],
        help="strain order (default: 2)",
    )
    parser.add_argument(
        "input_fn",
        type=str,
        nargs="?",
        const="POSCAR",
        default="POSCAR",
        help="input filename (default: POSCAR)",
    )
    parser.add_argument(
        "output_fn",
        type=str,
        nargs="?",
        const="strained.xyz",
        default="strained.xyz",
        help="output filename (default: strained.xyz)",
    )
    parser.add_argument(
        "-sl",
        "--strain_limit",
        type=float,
        default=0.1,
        help="strain limit",
    )

    args = parser.parse_args()

    strained_structure_generation(
        args.order,
        args.input_fn,
        args.output_fn,
        args.strain_limit,
    )
