"""采用 atomate 弛豫计算 workflow 默认参数"""

from atomate.common.powerups import add_namefile, add_tags
from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.workflows.presets.core import wf_structure_optimization
from fireworks.core.launchpad import LaunchPad
from pymatgen.core.structure import Structure

structure = Structure.from_prototype(
    prototype="diamond",
    species=["Si"],
    a=5.47,
)

user_incar_settings = {
    "ENCUT": 500,
    "ISPIN": 1,
    "ISMEAR": 0,
    "EDIFF": 1e-6,
    "EDIFFG": -1e-2,
}

wf = wf_structure_optimization(structure)

"""
# 解决 `cannot encode object: True, of type: <class 'numpy.bool'>` 报错
for key, val in wf.metadata.items():
    print(f"{key}: {type(val)}")

wf.metadata[key] = bool(wf.metadata[key])
"""

wf = add_modify_incar(wf, modify_incar_params={"incar_update": user_incar_settings})

wf = add_namefile(wf)
wf = add_tags(wf, {"task_name": "default relaxation workflow"})

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)

print("The default relaxation workflow is added.")
