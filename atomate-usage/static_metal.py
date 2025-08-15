"""atomate BCC Nb 静态计算（修改 INCAR、KPOINTS 参数以适用于金属体系）"""

from atomate.common.powerups import add_namefile, add_tags
from atomate.vasp.powerups import (
    add_modify_incar,
    add_modify_kpoints,
    add_modify_potcar,
)
from atomate.vasp.workflows.presets.core import wf_static
from fireworks.core.launchpad import LaunchPad
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints

structure = Structure.from_prototype(prototype="bcc", species=["Nb"], a=3.307)

user_incar_settings = (
    {
        "ISTART": 0,
        "ICHARG": 2,
        "ISPIN": 1,
        "LCHARG": True,
        "LWAVE": False,
        "PREC": "Accurate",
        "ENCUT": 500,
        "ALGO": "Normal",
        "ISMEAR": -5,
        "SIGMA": 0.05,
        "EDIFF": 1e-6,
        "NELM": 300,
        "NELMIN": 6,
        "NSW": 0,
        "IBRION": -1,
        "ISIF": 2,
    },
)

user_kpoints_settings = Kpoints(kpts=[[10, 10, 10]])

wf = wf_static(structure)

"""
# 解决 `cannot encode object: True, of type: <class 'numpy.bool'>` 报错
for key, val in wf.metadata.items():
    print(f"{key}: {type(val)}")

wf.metadata[key] = bool(wf.metadata[key])
"""

wf = add_modify_incar(wf, modify_incar_params={"incar_update": user_incar_settings})

wf = add_modify_kpoints(
    wf, modify_kpoints_params={"kpoints_update": user_kpoints_settings}
)

wf = add_modify_potcar(wf, modify_potcar_params={"potcar_symbols": {"Nb": "Nb_sv"}})

wf = add_namefile(wf)
wf = add_tags(wf, {"task_name": "Metal static workflow"})

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)

print("The Metal static workflow is added.")
