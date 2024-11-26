"""采用 atomate 弛豫计算 workflow 默认参数"""

from atomate.common.powerups import add_namefile, add_tags
from atomate.vasp.workflows.presets.core import wf_structure_optimization
from fireworks.core.launchpad import LaunchPad
from pymatgen.core.structure import Structure

structure = Structure.from_prototype(
    prototype="diamond",
    species=["Si"],
    a=5.47,
)
structure_copy = structure.copy()
structure_primitive = structure_copy.get_primitive_structure()

# structure_primitive = Structure.from_file("Si.vasp")

wf = wf_structure_optimization(structure_primitive)

"""
# 解决 `cannot encode object: True, of type: <class 'numpy.bool'>` 报错
for key, val in wf.metadata.items():
    print(f"{key}: {type(val)}")

wf.metadata[key] = bool(wf.metadata[key])
"""

wf = add_namefile(wf)
wf = add_tags(wf, {"task_name": "default relaxation workflow"})

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)

print("The default relaxation workflow is added.")
