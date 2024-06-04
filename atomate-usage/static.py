from atomate.common.powerups import add_namefile, add_tags
from atomate.vasp.workflows.presets.core import wf_static
from fireworks.core.launchpad import LaunchPad
from pymatgen.core.structure import Structure

# structure = Structure.from_prototype(prototype="diamond", species=["Si"], a=5.43)
# structure_copy = structure.copy()
# structure_primitive = structure_copy.get_primitive_structure()

structure = Structure.from_file("Si.vasp")

wf = wf_static(structure)

wf = add_namefile(wf)
wf = add_tags(wf, {"task_name": "atomate static workflow test"})

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)

print("The static test workflow is added.")
