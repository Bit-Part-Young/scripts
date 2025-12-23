#!/usr/bin/env python3

"""计算孪晶界面能"""

import argparse
import os

import numpy as np
from ase.io import read


def get_surface_area(structure_fn: str) -> float:
    """计算孪晶界面面积"""

    atoms = read(structure_fn, format="vasp")

    a, b, c = atoms.cell.lengths()
    alpha, beta, gamma = atoms.cell.angles()
    surface_area = a * b * np.sin(np.deg2rad(gamma))

    return surface_area


def get_twin_interface_energy2(twin_folder: str, gs_energy: float) -> float:
    """计算孪晶界面能（不构建并计算对应原子数的 bulk 构型）"""

    coeff = 1.6021766208 * 10**4

    structure_twin_fn = f"{twin_folder}/CONTCAR"
    surface_area_twin = get_surface_area(structure_twin_fn)

    atoms = read(f"{twin_folder}/OUTCAR", format="vasp-out")
    natoms = len(atoms)
    energy_twin = atoms.get_potential_energy()

    gamma = (energy_twin - natoms * gs_energy) / (2 * (surface_area_twin))
    gammaJ = gamma * coeff

    print(f"Twin interface energy: {gamma:.4f} eV/Å², {gammaJ:.1f} mJ/m².")


def get_twin_interface_energy(
    bulk_folder: str = "1-bulk", twin_folder: str = "2-twin"
) -> float:
    """计算孪晶界面能（构建并计算对应原子数的 bulk 构型）"""

    coeff = 1.6021766208 * 10**4

    structure_bulk_fn = f"{bulk_folder}/CONTCAR"
    structure_twin_fn = f"{twin_folder}/CONTCAR"
    if not os.path.exists(structure_bulk_fn):
        raise FileNotFoundError(f"File not found: {structure_bulk_fn}.")

    surface_area_bulk = get_surface_area(structure_bulk_fn)
    surface_area_twin = get_surface_area(structure_twin_fn)

    outcar_bulk = f"{bulk_folder}/OUTCAR"
    outcar_final = f"{twin_folder}/OUTCAR"

    if not os.path.exists(outcar_bulk):
        raise FileNotFoundError(f"File not found: {outcar_bulk}.")

    try:
        energy_initial = read(outcar_bulk, format="vasp-out").get_potential_energy()
        energy_final = read(outcar_final, format="vasp-out").get_potential_energy()
    except Exception as e:
        raise ValueError(f"Error reading OUTCAR file: {e}.")

    gamma = (energy_final - energy_initial) / (2 * (surface_area_bulk))
    gammaJ = gamma * coeff
    gamma2 = (energy_final - energy_initial) / (2 * (surface_area_twin))
    gammaJ2 = gamma2 * coeff

    print(
        f"\nEnergy: {bulk_folder} {energy_initial:.4f} eV, {twin_folder} {energy_final:.4f} eV."
    )
    print(
        f"Surface area: {bulk_folder} {surface_area_bulk:.4f} Å², {twin_folder} {surface_area_twin:.4f} Å²."
    )
    print(f"Twin interface energy (bulk area): {gamma:.4f} eV/Å², {gammaJ:.1f} mJ/m².")
    print(
        f"Twin interface energy (twin area): {gamma2:.4f} eV/Å², {gammaJ2:.1f} mJ/m²."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Twin interface energy calculation.")

    parser.add_argument("twin_folder", help="twin model calculation folder")
    parser.add_argument(
        "bulk",
        help="bulk model calculation folder or bulk unitcell ground state energy",
    )

    args = parser.parse_args()

    twin_folder = args.twin_folder
    bulk = args.bulk

    try:
        gs_energy = float(bulk)
        get_twin_interface_energy2(twin_folder, gs_energy)
    except:
        get_twin_interface_energy(bulk, twin_folder)
