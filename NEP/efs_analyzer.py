"""
reference: https://github.com/zhyan0603/GPUMDkit/blob/main/Scripts/analyzer/energy_force_virial_analyzer.py
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from ase.atoms import Atoms
from ase.io import read


def calculate_range(xyz_fn: str, property_name: str):
    """计算 energy, force, virial 的范围"""

    property_name = property_name.lower()
    # [ ] 待改进
    values = []

    atoms_list = read(xyz_fn, index=":", format="extxyz")

    for atoms in atoms_list:
        atoms: Atoms
        info_lower = {k.lower(): v for k, v in atoms.info.items()}

        if property_name == "energy":
            natoms = len(atoms)
            energy = atoms.get_potential_energy()
            energy_pa = energy / natoms
            values.append(energy_pa)
        elif property_name in ["force", "forces"]:
            forces = atoms.arrays["force"]
            values.extend(np.linalg.norm(forces, axis=1))
        elif property_name == "virial":
            if "virial" in info_lower:
                virial = info_lower["virial"]
                values.extend(virial)
            else:
                raise ValueError("Virial information not found in frame info.")
        else:
            raise ValueError("Invalid property. Choose from 'energy', 'force', or 'virial'.")

    return np.min(values), np.max(values), values


def plot_histogram(values, property_name: str):
    """绘制 force, energy, virial 的直方图"""

    property_name = property_name.capitalize()
    fig, ax = plt.subplots(figsize=(6, 4), dpi=300)

    ax.hist(values, bins=30, edgecolor="black")

    ax.set_title(f"{property_name} Histogram")
    ax.set_xlabel(f"{property_name}")
    ax.set_ylabel("Frequency")

    plt.tight_layout()

    fig.savefig(f"range_{property_name}.png")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate & plot the range of energy/forces/virials in an xyz file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "xyz_fn",
        type=str,
        nargs="?",
        default="train.xyz",
        const="train.xyz",
        help="xyz filename",
    )

    parser.add_argument(
        "-pn",
        "--property_name",
        type=str,
        choices=["energy", "forces", "virial"],
        help="Property name: energy, forces, or virial",
    )

    parser.add_argument(
        "--hist",
        action="store_true",
        help="Plot histogram of the property",
    )

    args = parser.parse_args()
    xyz_fn = args.xyz_fn
    property_name = args.property_name
    plot_hist = args.hist

    min_val, max_val, values = calculate_range(
        xyz_fn=xyz_fn,
        property_name=property_name,
    )

    print(f"{property_name.capitalize()} range: {min_val:.3f} to {max_val:.3f}")

    if plot_hist:
        plot_histogram(values, property_name)
