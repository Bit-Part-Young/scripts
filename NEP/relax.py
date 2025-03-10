"""结构弛豫 ASE 实现"""

import argparse

import numpy as np
from ase.atoms import Atoms
from ase.build import bulk
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS, FIRE, LBFGS, QuasiNewton


def setup_calculator(calculator_type: str):

    if calculator_type == "nep":
        from calorine.calculators import CPUNEP, GPUNEP

        return CPUNEP(model_filename="nep.txt")

    elif calculator_type == "mace-mpa":
        import torch
        from mace.calculators import mace_mp

        device = "cuda" if torch.cuda.is_available() else "cpu"

        return mace_mp(model="medium-mpa-0", device=device)

    elif calculator_type == "mattersim":
        import torch
        from mattersim.forcefield import MatterSimCalculator

        device = "cuda" if torch.cuda.is_available() else "cpu"

        return MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device=device)


def relax(
    atoms: Atoms,
    fmax: float = 0.01,
    steps: int = 500,
    model: str = "fire",
    constant_cell: bool = False,
    constant_volume: bool = False,
    **kwargs,
):
    if constant_cell:
        ecf = atoms
    else:
        ecf = FrechetCellFilter(atoms, constant_volume=constant_volume)

    if model == "qn":
        dyn = QuasiNewton(ecf, **kwargs)
    elif model == "bfgs":
        dyn = BFGS(ecf, **kwargs)
    elif model == "lbfgs":
        dyn = LBFGS(ecf, **kwargs)
    elif model == "fire":
        dyn = FIRE(ecf, **kwargs)

    dyn.run(fmax=fmax, steps=steps)

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    cell = atoms.cell[:]

    np.set_printoptions(precision=5, suppress=True)
    print(f"Energy: {round(energy, 6)}")
    print(f"Forces: {forces}")
    print(f"Cell: {cell}")


def main():
    atoms = bulk("Al", a=4.1, cubic=True)
    atoms = bulk("V", a=3.0, cubic=True)
    atoms = bulk("Ti", a=2.95, c=4.67)
    atoms = bulk("Zr", a=3.23, c=5.15)

    atoms.calc = setup_calculator(calculator_type="nep")

    relax(atoms)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Structure relaxation using ASE",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    main()
