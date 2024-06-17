from pymatgen.core.structure import Structure

"""
将不含 O 间隙原子的构型原子坐标及种类替换到含 O 原子的构型对应非 O 原子坐标及种类
"""


structure_tio = Structure.from_file("POSCAR_TiO")
structure_tialnb = Structure.from_file("POSCAR_TiAlNb")
num_oxygen = int(structure_tio.composition.as_dict()["O"])

structure_copy = structure_tio.copy()

for idx in range(structure_copy.num_sites - num_oxygen):
    structure_copy.replace(
        idx=idx,
        species=structure_tialnb[idx].specie,
        coords=structure_tialnb[idx].frac_coords,
    )

structure_copy.to(filename="POSCAR_TiO_replace")
