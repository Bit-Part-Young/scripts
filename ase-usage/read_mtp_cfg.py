from collections import deque

import numpy as np
from ase.parallel import paropen
from ase.utils import string2index


def get_max_index(index):
    if np.isscalar(index):
        return index
    elif isinstance(index, slice):
        return index.stop if (index.stop is not None) else float("inf")


# TODO: 待完成
def cfg_to_ase_atoms():
    """将 cfg 文件转换为 ase.Atoms 对象"""


# 写法参考 ase.io.lammpsrun 模块中的 ead_lammps_dump_text() 函数
def read_cfg(infileobj, index=-1):
    """读取 MTP 的 cfg 格式文件"""

    if isinstance(index, str):
        index = string2index(index)

    index_end = get_max_index(index)

    if isinstance(infileobj, str):
        fileobj = paropen(infileobj)

    lines = deque(fileobj.readlines())

    images = []

    natoms = None
    lattice = None
    array = None
    energy = None
    stress_list = None
    dft_code = None
    mindist = None
    while len(lines) > 0:
        line = lines.popleft()

        if "Size" in line:
            line = lines.popleft()
            natoms = int(line.strip())

        if "Supercell" in line:
            lattice = []
            for i in range(3):
                line = lines.popleft()
                lattice.append([float(x) for x in line.split()])
            lattice = np.array(lattice)

        # TODO: 待解析原子类型、坐标和受力
        if "AtomData" in line:
            array = []
            for i in range(natoms):
                line = lines.popleft()
                array.append([float(x) for x in line.split()])
            # print(positions)
            array = np.array(array)
            atoms_type = array[:, 1]
            positions_array = array[:, 2:5]
            force_array = array[:, 5:8]

        if "Energy" in line:
            line = lines.popleft()
            energy = float(line.strip())

        if "PlusStress" in line:
            line = lines.popleft()
            stress_list = [float(x) for x in line.split()]

        if "Feature   EFS_by" in line:
            dft_code = line.split()[2]
        if "Feature   mindist" in line:
            mindist = float(line.split()[2])
            data = {
                "natoms": natoms,
                "lattice": lattice,
                "atoms_type": atoms_type,
                "positions": positions_array,
                "forces": force_array,
                "energy": energy,
                "stress": stress_list,
                "dft_code": dft_code,
                "mindist": mindist,
            }

            images.append(data)

        if len(images) > index_end >= 0:
            break

    return images[index]


if __name__ == "__main__":

    # fn = "test.cfg"
    # fn = "test2.cfg"
    fn = "ZrSn.cfg"
    # read_cfg(fn)
    # test = read_cfg(infileobj=fn, index="8:10")
    test = read_cfg(infileobj=fn, index=-1)

    print(test)
    print(len(test))
