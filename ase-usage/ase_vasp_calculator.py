"""ase VASP Calculator 计算示例"""

from ase.build import molecule, bulk
from ase.calculators.vasp import Vasp

# N2 示例在 master 上运行会报错
# atoms = molecule("N2", pbc=True)
# atoms.center(vacuum=5)
atoms = bulk("Nb", a=3.2)


calc = Vasp(
    xc="pbe",  # 选择泛涵
    encut=400,  # 截断能
    # kpts=(1, 1, 1),  # k 点
    kpts=(5, 5, 5),  # k 点
)

# 执行计算
atoms.calc = calc
# 计算结束后，获取能量
energy = atoms.get_potential_energy()
print(f"Potential energy: {energy:.2f} eV.")
