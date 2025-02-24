"""结构弛豫 ASE 实现"""

import argparse

from ase.atoms import Atoms
from ase.build import bulk
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import FIRE, LBFGS, QuasiNewton


def setup_calculator(calculator_type: str, calculator_settings: dict):

    if calculator_type == "nep":
        from calorine.calculators import CPUNEP, GPUNEP

        return CPUNEP(model_filename="nep.txt")

    elif calculator_type == "mace-mpa":
        from mace.calculators import mace_mp

        return mace_mp(model="medium-mpa-0")

    elif calculator_type == "mattersim":
        from mattersim.forcefield import MatterSimCalculator

        return MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device="cpu")


def relax(atoms: Atoms, fmax: float = 0.01, steps: int = 500, model: str = "fire"):
    ecf = FrechetCellFilter(atoms)

    if model == "qn":
        dyn = QuasiNewton(ecf)
    elif model == "lbfgs":
        dyn = LBFGS(ecf)
    elif model == "fire":
        dyn = FIRE(ecf)

    dyn = FIRE(ecf)

    dyn.run(fmax=fmax, steps=steps)


def main():
    atoms = bulk("Mo", a=3.16, cubic=True)

    atoms.calc = setup_calculator()

    relax(atoms)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Structure relaxation using ASE",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    main()
