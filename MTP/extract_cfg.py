"""
提取 MTP cfg 文件 中的能量、力和应力数据

energy_train.out、force_train.out、stress_train.out 列的顺序对应 NEP 训练的输出文件
"""

import numpy as np
from ase.units import GPa


def extract_cfg(cfg_fn: str = "train.cfg") -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """提取 MTP cfg 文件 中的能量、力和应力数据"""

    energy_pa_list = []
    forces_list = []
    virial_list = []
    stress_list = []
    with open(cfg_fn) as f:
        line = "chongchongchong!"
        flag = 0
        while line:
            line = f.readline()
            if "BEGIN_CFG" in line:
                cell = np.zeros((3, 3))

            if "Size" in line:
                line = f.readline()
                natoms = int(line.split()[0])
                forces = np.zeros((natoms, 3))

            if "Supercell" in line:
                for i in range(3):
                    line = f.readline()
                    for j in range(3):
                        cell[i, j] = float(line.split()[j])

                volume = np.abs(np.linalg.det(cell))

            if "AtomData" in line:
                for i in range(natoms):
                    line = f.readline()
                    forces[i] = list(map(float, line.split()[-3:]))
                forces_list.append(forces)

            if "Energy" in line:
                line = f.readline()
                energy = float(line.split()[0])
                energy_pa = energy / natoms
                energy_pa_list.append(energy_pa)

            # PlusStress 实际对应 virial
            if "PlusStress" in line:
                line = f.readline()
                # MTP cfg stress 分量顺序 xx yy zz yz xz xy
                plus_stress_mtp = np.array(list(map(float, line.split())))
                # xyz stress 分量顺序 xx yy zz xy yz xz
                virial = [
                    plus_stress_mtp[0],
                    plus_stress_mtp[1],
                    plus_stress_mtp[2],
                    plus_stress_mtp[5],
                    plus_stress_mtp[3],
                    plus_stress_mtp[4],
                ]
                virial_pa = np.array(virial) / natoms
                virial_list.append(virial_pa)

                stress = -1 * (np.array(virial) / volume) / GPa
                stress_list.append(stress)

                flag += 1
                print(f"No. {flag} configuration processed!")

    energy_array = np.array(energy_pa_list).reshape(-1, 1)
    forces_array = np.concatenate(forces_list, axis=0)
    virial_array = np.array(virial_list)
    stress_array = np.array(stress_list)

    return energy_array, forces_array, virial_array, stress_array


def main():
    dft_cfg_fn = "train.cfg"
    mtp_cfg_fn = "train_predict.cfg"

    dft_energy, dft_forces, dft_virial, dft_stress = extract_cfg(dft_cfg_fn)
    mtp_energy, mtp_forces, mtp_virial, mtp_stress = extract_cfg(mtp_cfg_fn)

    energy_total = np.concatenate((mtp_energy, dft_energy), axis=1)
    forces_total = np.concatenate((mtp_forces, dft_forces), axis=1)
    virial_total = np.concatenate((mtp_virial, dft_virial), axis=1)
    stress_total = np.concatenate((mtp_stress, dft_stress), axis=1)

    energy_train_fn = "energy_train.out"
    force_train_fn = "force_train.out"
    virial_train_fn = "virial_total.out"
    stress_train_fn = "stress_train.out"
    np.savetxt(energy_train_fn, energy_total, fmt="%.10f")
    np.savetxt(force_train_fn, forces_total, fmt="%.10f")
    np.savetxt(virial_train_fn, virial_total, fmt="%.10f")
    np.savetxt(stress_train_fn, stress_total, fmt="%.10f")

    print("\nEnergy, Forces, Stress data shape:")

    print(
        energy_total.shape,
        forces_total.shape,
        virial_total.shape,
        stress_total.shape,
    )

    print("\nWork is done!")


if __name__ == "__main__":
    main()
