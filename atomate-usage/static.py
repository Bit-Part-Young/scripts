"""atomate Diamond Si 静态计算"""

from atomate.common.powerups import add_namefile, add_tags
from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.workflows.presets.core import wf_static
from fireworks.core.launchpad import LaunchPad
from pymatgen.core.structure import Structure

structure = Structure.from_prototype(prototype="diamond", species=["Si"], a=5.47)

user_incar_settings = {
    "ENCUT": 500,
    "ISPIN": 1,
    "ISMEAR": -5,
    "EDIFF": 1e-6,
}

wf = wf_static(structure)

"""
# 解决 `cannot encode object: True, of type: <class 'numpy.bool'>` 报错
for key, val in wf.metadata.items():
    print(f"{key}: {type(val)}")

wf.metadata[key] = bool(wf.metadata[key])
"""

wf = add_modify_incar(wf, modify_incar_params={"incar_update": user_incar_settings})

wf = add_namefile(wf)
wf = add_tags(wf, {"task_name": "default static workflow"})

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)

print("The default static workflow is added.")
