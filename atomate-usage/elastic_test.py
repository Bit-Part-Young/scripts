from atomate.common.powerups import add_namefile, add_tags
from atomate.vasp.powerups import add_modify_incar
from atomate.vasp.workflows.presets.core import wf_elastic_constant
from fireworks.core.launchpad import LaunchPad
from pymatgen.core.structure import Structure

structure = Structure.from_prototype(prototype="fcc", species=["Al"], a=4.401)

tags = {"tag_name": "elastic Al with atomate"}

user_incar_settings = {
    "ISMEAR": 1,
    "ISPIN": 1,
    "EDIFFG": -0.01,
}

wf = wf_elastic_constant(structure=structure)

wf = add_modify_incar(wf, modify_incar_params={"incar_update": user_incar_settings})

wf = add_namefile(wf)
wf = add_tags(wf, tags_list=tags)

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)

print("The elastic workflow is added.")
