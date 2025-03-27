"""
提取 MTP cfg 文件 中的能量、力和应力数据；不考虑位力

energy_train.out、force_train.out、stress_train.out 列的顺序对应 NEP 训练的输出文件
"""

import numpy as np


def extract_cfg(cfg_fn: str = "train.cfg") -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """提取 MTP cfg 文件 中的能量、力和应力数据"""

    energy_pa_list = []
    forces_list = []
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

            if "PlusStress" in line:
                line = f.readline()
                # MTP cfg stress 分量顺序 xx yy zz yz xz xy
                plusstress = np.array(list(map(float, line.split())))
                # xyz stress 分量顺序 xx yy zz xy yz xz
                nep_stress = [
                    plusstress[0],
                    plusstress[1],
                    plusstress[2],
                    plusstress[5],
                    plusstress[3],
                    plusstress[4],
                ]

                stress_list.append(nep_stress)

                flag += 1
                print(f"No. {flag} configuration processed!")

    energy_pa_array = np.array(energy_pa_list).reshape(-1, 1)
    forces_array = np.concatenate(forces_list, axis=0)
    stress_array = np.array(stress_list)

    return energy_pa_array, forces_array, stress_array


def main():
    dft_cfg_fn = "train.cfg"
    mtp_cfg_fn = "train_predict.cfg"

    dft_energy, dft_forces, dft_stress = extract_cfg(dft_cfg_fn)
    mtp_energy, mtp_forces, mtp_stress = extract_cfg(mtp_cfg_fn)

    # 预测值、DFT 计算值 列合并
    energy_concat = np.concatenate((mtp_energy, dft_energy), axis=1)
    forces_concat = np.concatenate((mtp_forces, dft_forces), axis=1)
    stress_concat = np.concatenate((mtp_stress, dft_stress), axis=1)

    energy_train_fn = "energy_train.out"
    force_train_fn = "force_train.out"
    stress_train_fn = "stress_train.out"
    np.savetxt(energy_train_fn, energy_concat, fmt="%.10f")
    np.savetxt(force_train_fn, forces_concat, fmt="%.10f")
    np.savetxt(stress_train_fn, stress_concat, fmt="%.10f")

    print("\nEnergy, Forces, Stress data shape:")

    print(
        energy_concat.shape,
        forces_concat.shape,
        stress_concat.shape,
    )

    print("\nWork is done!")


if __name__ == "__main__":
    main()
