"""将 json 构型及其数据文件转换为 ASE db 格式

文件格式: [{}, {}, ...]

dict_keys(['element', 'num_atoms', 'group', 'description', 'tag', 'structure', 'outputs'])

outputs 数据中的 key
dict_keys(['energy', 'forces', 'virial_stress'])
"""

import json

from ase.db import connect
from pymatgen.core.structure import Structure

db_fn = "quinary.db"
db = connect(db_fn)

# 文件共 9848 个构型
json_fn = "quinary.json"

with open(json_fn, "r") as f:
    json_data_list = json.load(f)


flag = 0
for json_data in json_data_list:
    structure = Structure.from_dict(json_data["structure"])
    atoms = structure.to_ase_atoms()

    del json_data["structure"]
    data = json_data["outputs"]
    del json_data["outputs"]
    key_value_pairs = json_data

    db.write(
        atoms,
        key_value_pairs=key_value_pairs,
        data=data,
    )

    flag += 1

    print(f"No. {flag} structure has been processed.")
