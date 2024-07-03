"""ASE VASP 能带结构计算"""

from ase.build import bulk
from ase.calculators.vasp import Vasp

si = bulk("Si")
# 计算目录
mydir = "bandstructure"

# 自洽计算
calc = Vasp(
    xc="PBE",
    kpts=(4, 4, 4),
    directory=mydir,
)

si.calc = calc
si.get_potential_energy()

print("calculation 1 is done!")

# 非自洽计算
kpts = {"path": "WGX", "npoints": 30}

calc.set(
    isym=0,  # Turn off kpoint symmetry reduction
    icharg=11,  # Non-SC calculation
    kpts=kpts,
)

si.get_potential_energy()

print("calculation 2 is done!")
print("\nwork is done!")
