#!/usr/bin/env python3

"""点缺陷形成能计算"""

import argparse

parser = argparse.ArgumentParser(description="Point defect formation energy", epilog="Author: SLY.")
parser.add_argument("pd_type", type=str, choices=["vac", "sia"], help="point defect type")
parser.add_argument("final_energy", type=float, help="Final energy with point defect")
parser.add_argument("bulk_energy", type=float, help="Bulk energy without point defect")
parser.add_argument("natoms", type=int, help="Number of atoms")

args = parser.parse_args()


def pde_cal(args):
    """点缺陷形成能计算"""

    pd_type = args.pd_type
    final_energy = args.final_energy
    bulk_energy = args.bulk_energy
    natoms = args.natoms

    if pd_type == "vac":
        point_defect_energy = final_energy - ((natoms - 1) / natoms) * bulk_energy
        print(
            f"Vacancy formation energy: {point_defect_energy:.2f} eV. Concentration: {1 / natoms * 100:.2f}%."
        )
    elif pd_type == "sia":
        point_defect_energy = final_energy - ((natoms + 1) / natoms) * bulk_energy
        print(
            f"SIA formation energy: {point_defect_energy:.2f} eV. Concentration: {1 / natoms * 100:.2f}%."
        )

    if point_defect_energy <= 0.0:
        print(
            "Warning: The calculated formation energy is negative, maybe the energy sequence is wrong. Please check!"
        )


pde_cal(args)

