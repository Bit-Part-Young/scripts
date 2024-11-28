"""采用 atomate 弛豫计算 workflow，修改 INCAR、KPOINTS 参数以适用于金属体系"""

from atomate.common.powerups import add_namefile, add_tags
from atomate.vasp.powerups import add_modify_incar, add_modify_kpoints
from atomate.vasp.workflows.presets.core import wf_structure_optimization
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints

from fireworks.core.launchpad import LaunchPad

structure = Structure.from_prototype(
    prototype="bcc",
    species=["Nb"],
    a=3.32,
)
structure_copy = structure.copy()
structure_primitive = structure_copy.get_primitive_structure()

# structure_primitive = Structure.from_file("Nb.vasp")

user_incar_settings = (
    {
        "ISTART": 0,
        "ICHARG": 2,
        "ISPIN": 1,
        "LCHARG": True,
        "LWAVE": False,
        "PREC": "Accurate",
        "ENCUT": 300,
        "ALGO": "Fast",
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "EDIFF": 1e-6,
        "NELM": 90,
        "NELMIN": 6,
        "NSW": 100,
        "IBRION": 1,
        "ISIF": 3,
        "EDIFFG": -1e-2,
    },
)

user_kpoints_settings = Kpoints.gamma_automatic(kpts=(5, 5, 5))

# 该方式 修改 INCAR 参数会报错
# c = {"USER_INCAR_SETTINGS": user_incar_settings}

wf = wf_structure_optimization(
    structure_primitive,
    # c=c,
)

"""
# 解决 `cannot encode object: True, of type: <class 'numpy.bool'>` 报错
for key, val in wf.metadata.items():
    print(f"{key}: {type(val)}")

wf.metadata[key] = bool(wf.metadata[key])
"""

wf = add_modify_incar(
    wf,
    modify_incar_params={
        "incar_update": user_incar_settings,
    },
)


wf = add_modify_kpoints(
    wf,
    modify_kpoints_params={
        "kpoints_update": user_kpoints_settings,
    },
)

wf = add_namefile(wf)
wf = add_tags(wf, {"task_name": "Metal relaxation workflow"})

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)

print("The Metal relaxation workflow is added.")
