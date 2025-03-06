"""
将 extxyz 文件格式转换为 cfg 文件格式

reference: https://github.com/hu-yanxiao/SUS2-MLIP/blob/main/python_tools/xyz2cfg.py
"""

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read


def set_spcalculator(atoms: Atoms):
    """设置 SinglePointCalculator 计算器"""

    energy = atoms.get_potential_energy()
    forces = atoms.arrays["force"]

    # xyz 格式中 virial 等同于 MTP cfg 中的 PlusStress
    # 用 virial 代替 stress
    virial = atoms.info["virial"]

    spcalculator = SinglePointCalculator(
        atoms=atoms,
        energy=energy,
        forces=forces,
        stress=virial,
    )

    atoms.calc = spcalculator

    return atoms


def write_cfg(
    xyz_fn: str,
    element_map: dict[str, int],
    cfg_fn: str = "output.cfg",
):
    """将 extxyz 文件格式转换为 cfg 文件格式"""

    atoms_list = read(xyz_fn, index=":", format="extxyz")

    ff = open(cfg_fn, "w")
    for atoms in atoms_list:
        atoms = set_spcalculator(atoms)

        chemical_symbols = atoms.get_chemical_symbols()
        natoms = len(chemical_symbols)

        cell = atoms.get_cell()
        positioins = atoms.get_positions()

        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        plusstress = atoms.get_stress(voigt=False)

        ff.write("""BEGIN_CFG\n""")
        ff.write(""" Size\n""")
        ff.write("""  {:6}\n""".format(natoms))
        ff.write(""" Supercell \n""")
        ff.write("""{:15.10f} {:15.10f} {:15.10f}\n""".format(cell[0, 0], cell[0, 1], cell[0, 2]))
        ff.write("""{:15.10f} {:15.10f} {:15.10f}\n""".format(cell[1, 0], cell[1, 1], cell[1, 2]))
        ff.write("""{:15.10f} {:15.10f} {:15.10f}\n""".format(cell[2, 0], cell[2, 1], cell[2, 2]))
        ff.write(
            """ AtomData:  id type       cartes_x      cartes_y      cartes_z     fx          fy          fz\n"""
        )
        for atom_index in range(natoms):
            ff.write(
                """ {:6} {:6} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f}\n""".format(
                    atom_index + 1,
                    element_map[chemical_symbols[atom_index]],
                    positioins[atom_index, 0],
                    positioins[atom_index, 1],
                    positioins[atom_index, 2],
                    forces[atom_index, 0],
                    forces[atom_index, 1],
                    forces[atom_index, 2],
                )
            )
        ff.write(""" Energy \n""")
        ff.write(f"""     {energy:15.10f} \n""")
        ff.write(
            """ PlusStress:  xx          yy          zz          yz          xz          xy\n"""
        )
        ff.write(
            "     {:15.10f} {:15.10f} {:15.10f} {:15.10} {:15.10f} {:15.10}\n".format(
                plusstress[0, 0],
                plusstress[1, 1],
                plusstress[2, 2],
                plusstress[1, 2],
                plusstress[0, 2],
                plusstress[0, 1],
            )
        )
        ff.write("""END_CFG\n""")
    ff.close()


if __name__ == "__main__":
    xyz_fn = "test.xyz"
    element_map = {
        "Al": 0,
        "Ti": 1,
        "V": 2,
        "Zr": 3,
    }
    write_cfg(xyz_fn, element_map=element_map, cfg_fn="output.cfg")
