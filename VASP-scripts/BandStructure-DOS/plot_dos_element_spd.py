"""元素分态密度绘制（使用 for 循环）"""

import matplotlib.pyplot as plt
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.outputs import Vasprun
from spt.plot_params import set_roman_plot_params

dos_vasprun = Vasprun("./dos/vasprun.xml")
dos_data = dos_vasprun.complete_dos

# 存储元素分态密度
element_pdos_dict = {}
element_orbitals = {
    "Nb": [
        OrbitalType.s,
        OrbitalType.p,
        OrbitalType.d,
    ],
    "Si": [
        OrbitalType.s,
        OrbitalType.p,
        # OrbitalType.d,
    ],
}

for element, orbitals in element_orbitals.items():
    element_pdos_dict[element] = {}
    for orbital in orbitals:
        element_pdos_dict[element][str(orbital)] = (
            dos_data.get_element_spd_dos(element)[orbital]
        )
# print(pdos_dict)


# 态密度绘制
set_roman_plot_params()

dos_plotter = DosPlotter()
# 总态密度
dos_plotter.add_dos("Total", dos_data)
# 元素分态密度
for element, orbitals in element_pdos_dict.items():
    for orbital, pdos in orbitals.items():
        pdos_data = dos_plotter.add_dos(f"{element}({orbital})", pdos)

ax = dos_plotter.get_plot(
    xlim=(-10, 2.5),
    ylim=(0, 35),
)

plt.savefig("dos_element_spd.png")
