"""根据构型生成 K-Path"""

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath

structure = Structure.from_file("POSCAR")

kpath = HighSymmKpath(structure)

# division = 20 vaspkit 默认值
kpoints = Kpoints.automatic_linemode(divisions=20, ibz=kpath)

print(kpoints)
