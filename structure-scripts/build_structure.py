"""
构建晶体结构

reference: [python学习之创建LAMMPS可读取的晶体结构模型](https://mp.weixin.qq.com/s/6bkomDigE4krXxc80-Oaiw)
"""

import numpy as np

# FCC Al 晶格信息
lattice_constant = 4.05
basis_vector = (
    np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    * lattice_constant
)
base_atoms = (
    np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ]
    )
    * lattice_constant
)

# 扩胞
# 超胞尺寸
system_size = 2

# 写入超胞原子坐标
positions = []
for i in range(system_size):
    for j in range(system_size):
        for k in range(system_size):
            base_position = np.array([i, j, k])
            cartesian_position = np.inner(basis_vector.T, base_position)
            print(cartesian_position)
            for atom in base_atoms:
                positions.append(cartesian_position + atom)

# 写入 LAMMPS data 文件
with open("data.lmp", "w") as fdata:
    fdata.write("Crystalline Al atoms - written for EnCodeVentor tutorial\n\n")
    fdata.write("{} atoms\n".format(len(positions)))
    fdata.write("{} atom types\n".format(1))
    fdata.write("{} {} xlo xhi\n".format(0.0, system_size * lattice_constant))
    fdata.write("{} {} ylo yhi\n".format(0.0, system_size * lattice_constant))
    fdata.write("{} {} zlo zhi\n".format(0.0, system_size * lattice_constant))
    fdata.write("\n")
    fdata.write("Atoms\n\n")
    for i, pos in enumerate(positions):
        fdata.write("{} 1 {} {} {}\n".format(i + 1, *pos))
