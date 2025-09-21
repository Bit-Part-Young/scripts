#!/usr/bin/env python3

"""使用 dpdata 生成随机 cell & position 扰动结构"""

import argparse

import dpdata
from ase.io import write


def perturb_dpdata(
    input_fn: str = "POSCAR",
    pert_num: int = 50,
    cell_pert_fraction: float = 0.10,
    atom_pert_distance: float = 0.1,
    output_fn: str = "perturbed.xyz",
):
    """使用 dpdata 生成随机 cell & position 扰动结构"""

    perturbed_system = dpdata.System(input_fn, fmt="vasp/poscar").perturb(
        pert_num=pert_num,
        cell_pert_fraction=cell_pert_fraction,
        atom_pert_distance=atom_pert_distance,
        atom_pert_style="uniform",
        atom_pert_prob=1.0,
    )

    atoms_list = perturbed_system.to("ase/structure")
    write(output_fn, atoms_list, format="extxyz", append=True)

    print(f"Total {pert_num} perturbed structures written to {output_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate random cell & position perturbed structures using dpdata.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "input_fn",
        type=str,
        default="POSCAR",
        help="input structure filename (default: POSCAR)",
    )
    parser.add_argument(
        "output_fn",
        type=str,
        default="perturbed.xyz",
        help="output xyz filename (default: perturbed.xyz)",
    )
    parser.add_argument(
        "-pn",
        "--pert_num",
        type=int,
        default=50,
        metavar="N",
        help="number of perturbed structures",
    )
    parser.add_argument(
        "-cp",
        "--cell_perturb",
        type=float,
        default=0.10,
        metavar="FLOAT",
        help="cell perturbation fraction",
    )
    parser.add_argument(
        "-ap",
        "--atom_perturb",
        type=float,
        default=0.1,
        metavar="FLOAT",
        help="atom perturbation distance magnitude",
    )

    args = parser.parse_args()

    perturb_dpdata(
        input_fn=args.input_fn,
        output_fn=args.output_fn,
        pert_num=args.pert_num,
        cell_pert_fraction=args.cell_perturb,
        atom_pert_distance=args.atom_perturb,
    )
