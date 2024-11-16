"""获取体系分态密度数据"""

import pandas as pd
from pymatgen.electronic_structure.core import OrbitalType, Spin
from pymatgen.io.vasp.outputs import Vasprun

dos_vasprun = Vasprun("./dos/vasprun.xml")
dos_data = dos_vasprun.complete_dos

# 获取费米能级
fermi = dos_data.efermi
# 整体能量平移
energy = dos_data.energies - fermi

# 存储体系分态密度
pdos_densities_dict = {}
orbital_list = [
    OrbitalType.s,
    OrbitalType.p,
    OrbitalType.d,
]

spd_dos = dos_data.get_spd_dos()
for orbital in orbital_list:
    dos_i = spd_dos[orbital]
    # 体系分态密度数据
    densities = dos_i.densities[Spin.up]
    pdos_densities_dict[str(orbital)] = densities

energy_dos_dict = {
    "energy": energy,
    "Total": dos_data.densities[Spin.up],
    **pdos_densities_dict,
}


df = pd.DataFrame(energy_dos_dict)
# print(df)
print(f'Max density: {df["energy"].max()}')

csv_fname = "dos_spd.csv"
df.to_csv(csv_fname, index=False)
