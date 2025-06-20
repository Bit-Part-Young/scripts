"""统计多晶模型中晶界和 Bulk 区域的元素分布及占比"""

import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

import numpy as np
import pandas as pd
from ovito.io import import_file
from ovito.modifiers import (
    CommonNeighborAnalysisModifier,
    DeleteSelectedModifier,
    SelectTypeModifier,
)

# pipeline = import_file("2-MCMD/restart.xyz")
# pipeline = import_file("2-MCMD/dump.xyz")
# pipeline = import_file("2-MCMD/backup/model.xyz")
# pipeline = import_file("2-MCMD/backup/restart.xyz")
# pipeline = import_file("2-MCMD/backup/dump.xyz")
# pipeline = import_file("1-Release-Stress/restart.xyz")
pipeline = import_file("1-Release-Stress/model.xyz")
print(pipeline.num_frames)

cna_modifier = CommonNeighborAnalysisModifier()
pipeline.modifiers.append(cna_modifier)

select_type_modifier = SelectTypeModifier(
    operate_on="particles",
    property="Structure Type",
    types={
        # Bulk 区域
        # CommonNeighborAnalysisModifier.Type.OTHER,
        # 晶界区域
        CommonNeighborAnalysisModifier.Type.BCC,
        CommonNeighborAnalysisModifier.Type.FCC,
        CommonNeighborAnalysisModifier.Type.HCP,
        CommonNeighborAnalysisModifier.Type.ICO,
    },
)
pipeline.modifiers.append(select_type_modifier)

pipeline.modifiers.append(DeleteSelectedModifier())

print(pipeline.modifiers)

element_info_list = []
for frame in range(pipeline.num_frames):
    # if frame < 600:
    #     continue
    data = pipeline.compute(frame)

    print(f"No. {frame} frame processed.")

    partcle_type_array = data.particles["Particle Type"]
    type_count_tuple = np.unique(partcle_type_array, return_counts=True)
    # {1: 100, ...}
    type_count_dict = dict(zip(type_count_tuple[0], type_count_tuple[1]))

    # {"frame": 1, "Ti": 100, ...}
    element_count_dict = {"frame": frame}
    for type in data.particles.particle_types.types:
        element_count_dict[type.name] = type_count_dict[type.id]

    element_info_list.append(element_count_dict)

    # bcc_count = data.attributes["CommonNeighborAnalysis.counts.BCC"]
    # fcc_count = data.attributes["CommonNeighborAnalysis.counts.FCC"]
    # hcp_count = data.attributes["CommonNeighborAnalysis.counts.HCP"]
    # ico_count = data.attributes["CommonNeighborAnalysis.counts.ICO"]
    # other_count = data.attributes["CommonNeighborAnalysis.counts.OTHER"]

    # print(frame, bcc_count, fcc_count, hcp_count, ico_count, other_count)


df = pd.DataFrame(element_info_list)
print(df)

df_array = df.to_numpy()
element_count_mean = df_array[:, 1:].mean(axis=0).astype(int)
print(element_count_mean)
element_count_mean_ratio = (element_count_mean / element_count_mean.sum()).round(2)
print(element_count_mean_ratio)

element_count_std = df_array[:, 1:].std(axis=0).astype(int)
print(element_count_std)

# df.to_csv("gb_particle_info.csv", index=False, sep=" ")
# df.to_csv("bulk_particle_info.csv", index=False, sep=" ")
