"""使用 pymatgen 模块进行能带绘制"""

from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp.outputs import Vasprun
from spt.plot_params import set_roman_plot_params

bs_vasprun = Vasprun("vasprun.xml", parse_projected_eigen=True)
bs_data = bs_vasprun.get_band_structure(line_mode=True)

set_roman_plot_params()

bs_plot = BSPlotter(bs=bs_data)
ax = bs_plot.get_plot(ylim=(-3, 3))
fig = ax.get_figure()
# 调整能带图的大小
fig.set_size_inches(6, 10)

ax.set_xticklabels(["M", "$\Gamma$", "K", "M"])
ax.legend([])

fig.savefig("bs_pymatgen.jpg")

print("Band structure figure is generated.")
