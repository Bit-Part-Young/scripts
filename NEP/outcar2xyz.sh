import dpdata
import numpy as np
from calorine.tools import relax_structure

d_vasprun = dpdata.LabeledSystem("vasprun.xml")

print(d_vasprun)
print(d_vasprun["energies"])
print(d_vasprun["forces"])
print(d_vasprun["virials"])
# print(d_vasprun["stress"])
# print(d_vasprun["cells"])

for cell in d_vasprun["cells"]:
    volume = np.abs(np.linalg.det(cell))
    print(volume)

supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
