"""ase VASP Caculator 计算示例"""

from ase.build import molecule
from ase.calculators.vasp import Vasp

atoms = molecule("N2", pbc=True)
atoms.center(vacuum=5)


calc = Vasp(
    xc="pbe",  # 选择泛涵
    encut=400,  # 截断能
    kpts=(1, 1, 1),  # k 点
)

# 执行计算
atoms.calc = calc
# 计算结束后，获取能量
energy = atoms.get_potential_energy()
print(f"Potential energy: {energy:.2f} eV.")
