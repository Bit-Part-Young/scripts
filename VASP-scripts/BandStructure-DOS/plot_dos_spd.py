"""体系分态密度绘制（使用 for 循环）"""

import matplotlib.pyplot as plt
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.outputs import Vasprun
from spt.plot_params import set_roman_plot_params

dos_vasprun = Vasprun("./dos/vasprun.xml")
dos_data = dos_vasprun.complete_dos

# 存储体系分态密度
pdos_dict = {}
orbital_list = [
    OrbitalType.s,
    OrbitalType.p,
    OrbitalType.d,
]

spd_dos = dos_data.get_spd_dos()
for orbital in orbital_list:
    pdos_dict[str(orbital)] = spd_dos[orbital]
# print(pdos_dict)


# 态密度绘制
set_roman_plot_params()

dos_plotter = DosPlotter()
# 总态密度
dos_plotter.add_dos("Total", dos_data)
# 分态密度
for orbital, pdos in pdos_dict.items():
    pdos_data = dos_plotter.add_dos(f"{orbital}", pdos)

ax = dos_plotter.get_plot(
    xlim=(-10, 2.5),
    ylim=(0, 35),
)

plt.savefig("dos_spd.png")

print("DOS Figure is generated.")
