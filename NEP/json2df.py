"""
将 json 构型及其数据文件转换为 Pandas DataFrame 格式

文件格式: [{}, {}, ...]

dict_keys(['element', 'num_atoms', 'group', 'description', 'tag', 'structure', 'outputs'])

outputs 数据中的 key
dict_keys(['energy', 'forces', 'virial_stress'])
"""

import pandas as pd
from monty.serialization import loadfn

json_fn = "quinary.json"
json_data_list = loadfn(json_fn)

general_info_list = []

flag = 0
for json_data in json_data_list:

    # 将 structure 和 outputs 中的 force virial_stress 数据删除
    # 其不易保存成 csv 文件
    del json_data["structure"]
    energy = json_data["outputs"]["energy"]
    del json_data["outputs"]
    general_info_dict = json_data
    general_info_dict["energy"] = energy
    general_info_list.append(general_info_dict)

    flag += 1

    print(f"No. {flag} structure has been processed.")

df = pd.DataFrame(general_info_list)
csv_fn = "quinary_info.csv"
df.to_csv(csv_fn, index=False)
