"""
读取 MTP 的 cfg 格式文件并转换为 ase.Atoms 对象

Author: SLY
Date: 2024-12-08
Updated: 2025-03-27
"""

from collections import deque

import numpy as np
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.parallel import paropen
from ase.units import GPa
from ase.utils import string2index


def get_max_index(index):
    if np.isscalar(index):
        return index
    elif isinstance(index, slice):
        return index.stop if (index.stop is not None) else float("inf")


def cfg2Atoms(
    symbols,
    cell,
    positions,
    pbc=True,
    info=None | dict,
    forces=None | np.ndarray,
    energy=None | float,
):
    """将 cfg 文件转换为 ase.Atoms 对象"""

    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=pbc,
        info=info,
    )

    atoms.arrays["forces"] = forces

    volume = atoms.get_volume()
    # -1 是为了与 ASE 做法一致
    virial = np.array(info["virial"])
    # stress 需是 6 个分量且为 voigt order
    stress = -1 * virial / volume

    # 将 virial 重写成 9 个分量
    virial = np.array(
        [
            [virial[0], virial[5], virial[4]],
            [virial[5], virial[1], virial[3]],
            [virial[4], virial[3], virial[2]],
        ]
    )
    atoms.info["virial"] = virial

    if forces is not None:
        calculator = SinglePointCalculator(
            atoms=atoms,
            energy=energy,
            forces=forces,
            stress=stress,
        )
        atoms.calc = calculator

    return atoms


# 写法参考 ase.io.lammpsrun 模块中的 read_lammps_dump_text() 函数
def read_mtp_cfg(
    infileobj: str,
    index=-1,
    symbols_map=dict[int, str],
):
    """读取 MTP 的 cfg 格式文件"""

    if isinstance(index, str):
        index = string2index(index)

    index_end = get_max_index(index)

    if isinstance(infileobj, str):
        fileobj = paropen(infileobj)

    lines = deque(fileobj.readlines())

    images = []

    while len(lines) > 0:
        line = lines.popleft()

        if "Size" in line:
            line = lines.popleft()
            natoms = int(line.strip())

        if "Supercell" in line:
            cell = []
            for _ in range(3):
                line = lines.popleft()
                cell.append([float(x) for x in line.split()])
            cell = np.array(cell)

        if "AtomData" in line:
            atomdata_list = []
            for _ in range(natoms):
                line = lines.popleft()
                atomdata_list.append([float(x) for x in line.split()])

            atomdata_array = np.array(atomdata_list)

            symbols_int = atomdata_array[:, 1].tolist()
            symbols = [symbols_map[symbol_int] for symbol_int in symbols_int]

            positions = atomdata_array[:, 2:5]
            forces = atomdata_array[:, 5:8]

        if "Energy" in line:
            line = lines.popleft()
            energy = float(line.strip())

        if "PlusStress" in line:
            line = lines.popleft()
            # MTP 用 VASP OUTCAR 中的 Total 作为其 cfg 构型文件中的 PlusStress 数据
            # 6 个分量时，MTP 和 ASE 采用的都是 voigt order
            plusstress = list(map(float, line.split()))
            virial_dict = {"virial": plusstress}

            atoms = cfg2Atoms(
                symbols=symbols,
                cell=cell,
                positions=positions,
                energy=energy,
                forces=forces,
                info=virial_dict,
            )
            images.append(atoms)

        if "Feature   EFS_by" in line:
            dft_code = line.split()[2]
            atoms.info["dft_code"] = dft_code

        # 该 Feature 不一定有
        if "Feature   mindist" in line:
            mindist = float(line.split()[2])
            atoms.info["mindist"] = mindist

        if len(images) > index_end >= 0:
            break

    return images[index]


if __name__ == "__main__":
    cfg_fn = "example.cfg"

    symbols_map = {0: "Zr", 1: "Sn"}

    atoms_list = read_mtp_cfg(cfg_fn, index=":", symbols_map=symbols_map)

    print(f"Total {len(atoms_list)} configurations.\n")

    atoms: Atoms = atoms_list[0]

    np.set_printoptions(precision=5, suppress=True)
    print(f"No. 1 configuration:\n {atoms}")
    print(f"Info:\n{atoms.info}")
    print(f"Energy (eV): {atoms.get_potential_energy()}")
    print(f"Forces (eV/Å):\n{atoms.get_forces()[:3]}")
    print(f"Stress (eV/Å^3):\n{atoms.get_stress(voigt=False)}")
    print(f"Stress (GPa):\n{atoms.get_stress(voigt=False) / GPa}")
    print(f"Virial (eV):\n{atoms.info['virial']}")
