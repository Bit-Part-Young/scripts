"""
读取 MTP 的 cfg 格式文件并转换为 ase.Atoms 对象

Author: SLY
Date: 2024-12-08
Updated: 2025-02-15
"""

from collections import deque

import numpy as np
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.parallel import paropen
from ase.utils import string2index


def get_max_index(index):
    if np.isscalar(index):
        return index
    elif isinstance(index, slice):
        return index.stop if (index.stop is not None) else float("inf")


def cfg_to_atoms(
    symbols,
    cell,
    positions,
    info=None,
    pbc=True,
    forces=None,
    energy=None,
):
    """将 cfg 文件转换为 ase.Atoms 对象"""

    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=pbc,
        info=info,
    )

    volume = atoms.get_volume()
    stress = -np.array(info["virial"]) / volume

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

    natoms = None
    cell = None
    array_total = None
    energy = None
    plusstress = None
    dft_code = None
    mindist = None
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
            array_total = []
            for _ in range(natoms):
                line = lines.popleft()
                array_total.append([float(x) for x in line.split()])
            array_total = np.array(array_total)
            atoms_type = array_total[:, 1]
            positions = array_total[:, 2:5]
            forces = array_total[:, 5:8]

        if "Energy" in line:
            line = lines.popleft()
            energy = float(line.strip())

        if "PlusStress" in line:
            line = lines.popleft()
            plusstress = [float(x) for x in line.split()]
            virial = [
                plusstress[0],
                plusstress[5],
                plusstress[4],
                plusstress[5],
                plusstress[1],
                plusstress[3],
                plusstress[4],
                plusstress[3],
                plusstress[2],
            ]
            virial_dict = {"virial": " ".join(map(str, virial))}

            symbols = [symbols_map[x] for x in atoms_type]
            out_atoms = cfg_to_atoms(
                symbols=symbols,
                cell=cell,
                positions=positions,
                energy=energy,
                forces=forces,
                info=virial_dict,
            )
            images.append(out_atoms)

        if "Feature   EFS_by" in line:
            dft_code = line.split()[2]

            out_atoms.info["dft_code"] = dft_code

        # 该 Feature 不一定有
        if "Feature   mindist" in line:
            mindist = float(line.split()[2])

            out_atoms.info["mindist"] = mindist

        if len(images) > index_end >= 0:
            break

    return images[index]


if __name__ == "__main__":
    cfg_fn = "example.cfg"

    symbols_map = {0: "Zr", 1: "Sn"}

    atoms_list = read_mtp_cfg(cfg_fn, index=":", symbols_map=symbols_map)

    print(f"Total {len(atoms_list)} configurations.\n")

    atoms: Atoms = atoms_list[0]

    print(f"No. 1 configuration:\n {atoms}")
    print(f"Info: {atoms.info}")
    print(f"Energy: {atoms.get_potential_energy()}")
    print(f"Forces:\n{atoms.get_forces()}")
    print(f"Stress:\n{atoms.get_stress()}")
    print(f"Virial:\n{atoms.info['virial']}")
