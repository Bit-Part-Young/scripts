"""
修改 atomate 源代码
"""

from atomate.vasp.config import (
    ADD_WF_METADATA,
    DB_FILE,
    SMALLGAP_KPOINT_MULTIPLY,
    STABILITY_CHECK,
    VASP_CMD,
)
from atomate.vasp.powerups import (
    add_common_powerups,
    add_modify_incar,
    add_small_gap_multiply,
    add_stability_check,
    add_wf_metadata,
)
from atomate.vasp.workflows.base.core import get_wf
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPHSERelaxSet, MPRelaxSet, MPStaticSet

# --------------------------------------------------------------


# atomate/vasp/workflows/presets/core.py
def wf_static(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    user_incar_settings = c.get("USER_INCAR_SETTINGS")
    user_kpoints_settings = c.get("USER_KPOINTS_SETTINGS")

    wf = get_wf(
        structure,
        "static_only.yaml",
        vis=MPStaticSet(
            structure,
            force_gamma=True,
            user_incar_settings=user_incar_settings,
            user_kpoints_settings=user_kpoints_settings,
        ),
        common_params={"vasp_cmd": vasp_cmd, "db_file": db_file},
    )

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_structure_optimization(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    user_incar_settings = c.get("USER_INCAR_SETTINGS")
    user_kpoints_settings = c.get("USER_KPOINTS_SETTINGS")

    wf = get_wf(
        structure,
        "optimize_only.yaml",
        vis=MPRelaxSet(
            structure,
            force_gamma=True,
            user_incar_settings=user_incar_settings,
            user_kpoints_settings=user_kpoints_settings,
        ),
        common_params={"vasp_cmd": vasp_cmd, "db_file": db_file},
    )

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


# --------------------------------------------------------------
