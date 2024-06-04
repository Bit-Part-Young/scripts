from atomate.common.powerups import add_namefile, add_tags
from atomate.vasp.workflows.presets.core import wf_structure_optimization
from fireworks.core.launchpad import LaunchPad
from pymatgen.core.structure import Structure

structure = Structure.from_file("Si.vasp")

wf = wf_structure_optimization(structure)

wf = add_namefile(wf)
wf = add_tags(wf, {"task_name": "atomate relaxation workflow test"})

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)

print("The relaxation test workflow is added.")
