"""
Phonopy 声子谱计算的 Python 脚本

reference: https://github.com/huiju-lee/HT-Phonon-MLIP/blob/main/calc_phonon.py
"""

import phonopy
from ase.calculators import calculator
from ase.calculators.vasp import Vasp
from ase.dft.kpoints import get_special_points
from mace.calculators import MACECalculator, mace_mp
from phonopy import Phonopy
from phonopy.file_IO import (
    parse_disp_yaml,
    parse_FORCE_CONSTANTS,
    parse_FORCE_SETS,
    write_disp_yaml,
    write_disp_yaml_from_dataset,
    write_FORCE_CONSTANTS,
    write_FORCE_SETS,
)
from phonopy.interface.calculator import read_crystal_structure
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core.structure import Structure

# 读取结构，转换成 PhonopyAtoms 格式的 unitcell
structure_fn = "POSCAR"
unitcell, _ = read_crystal_structure(structure_fn, interface_mode="vasp")

print(unitcell)
"""
lattice:
- [     3.320000000000000,     0.000000000000000,     0.000000000000000 ] # a
- [     0.000000000000000,     3.320000000000000,     0.000000000000000 ] # b
- [     0.000000000000000,     0.000000000000000,     3.320000000000000 ] # c
points:
- symbol: Nb # 1
  coordinates: [  0.000000000000000,  0.000000000000000,  0.000000000000000 ]
  mass: 92.906380
- symbol: Nb # 2
  coordinates: [  0.500000000000000,  0.500000000000000,  0.500000000000000 ]
  mass: 92.906380
"""

# 初始化 Phonopy
supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
phonon = Phonopy(
    unitcell=unitcell,
    supercell_matrix=supercell_matrix,
)

# 生成含位移的超胞
phonon.generate_displacements(distance=0.03)  # 默认值 0.01
swds = phonon.supercells_with_displacements
print(swds)  # phonopy.structure.atoms.PhonopyAtoms 列表
print(len(swds))
print(swds[0].cell, swds[0].scaled_positions, swds[0].symbols)
dataset = phonon.get_displacement_dataset()

displacement_data = []
for swd in swds:
    structure = Structure(
        lattice=swd.cell,
        species=swd.symbols[:],
        coords=swd.scaled_positions,
    )
    atoms = structure.to_ase_atoms()
    MACECalculator.calculate(atoms=atoms, properties="forces")
    forces = calculator.results["forces"]
    displacement_data.append({"structure": structure, "forces": forces})


def write_force_constants():
    phonon = phonopy.load()
    phonon.produce_force_constants()
    phonon.force_constants
    phonon.save()


def misc():
    phonon = phonopy.load()
    phonon.auto_band_structure()


def write_thermal_properties_yaml():
    phonon = phonopy.load()
    phonon.run_mesh()
    phonon.run_thermal_properties()
    phonon.write_yaml_thermal_properties()


def misc():
    phonon = phonopy.load()
    phonon.auto_band_structure()
    phonon.run_total_dos()
    phonon.run_projected_dos()
    phonon.get_band_structure_dict()
    phonon.get_band_structure()
    phonon.get_total_DOS()
    phonon.get_partial_DOS()
    phonon.thermal_properties.write_yaml()
