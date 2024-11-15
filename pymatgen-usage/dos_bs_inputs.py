# 生成 DOS 和能带结构计算的输入文件

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPNonSCFSet, MPRelaxSet, MPStaticSet
from pymatgen.symmetry.bandstructure import HighSymmKpath

structure_fn = ...
structure = Structure.from_file(structure_fn)

# 结构优化
incar_relax = {
    "ENCUT": 400,
    "EDIFF": 1e-4,
    "EDIFG": -0.05,
    "ISPIN": 1,
    "LORBIT": 12,
}

relax_settings = MPRelaxSet(
    structure,
    user_incar_settings=incar_relax,
    user_kpoints_settings=Kpoints.gamma_automatic(kpts=...),
)
relax_settings.write_input("./relax")
print("Relax calculation input files is generated!")

# 静态/自洽计算
static_settings = MPStaticSet.from_prev_calc(
    prev_calc_dir="./relax",
)
static_settings.write_input("./static")
print("Static calculation input files is generated!")

# 非自洽计算
nonscf_dos = MPNonSCFSet.from_prev_calc(
    prev_calc_dir="./static",
    user_kpoints_settings=Kpoints.gamma_automatic(kpts=...),
)
nonscf_dos.write_input("./dos")
print("DOS calculation input files is generated!")

# 能带结构计算
# 生成高对称路径
structure_prim = structure.get_primitive_structure()
kpath = HighSymmKpath(structure=structure_prim, path_type="hinuma")
# 根据高对称路径生成 K 点
kpoints_bs = Kpoints.automatic_linemode(divisions=10, ibz=kpath)
incar_bs = {
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "LORBIT": 11,
    "NBANDS": 64,
}

bs_settings = MPNonSCFSet.from_prev_calc(
    prev_calc_dir="./static",
    user_incar_settings=incar_bs,
    user_kpoints_settings=kpoints_bs,
)
bs_settings.write_input("./bs")

print("DOS calculation input files is generated!")
