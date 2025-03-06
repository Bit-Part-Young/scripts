"""提取 MTP cfg 文件 中的能量、力和应力数据"""

import numpy as np


def extract_cfg(cfg_fn: str = "train.cfg") -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """提取 MTP cfg 文件 中的能量、力和应力数据"""

    forces_list = []
    energy_list = []
    stress_list = []
    with open(cfg_fn) as f:
        line = "chongchongchong!"
        flag = 0
        while line:
            line = f.readline()
            if "Size" in line:
                line = f.readline()
                natoms = int(line.split()[0])
                forces = np.zeros((natoms, 3))

            if "AtomData" in line:
                for i in range(natoms):
                    line = f.readline()
                    forces[i] = list(map(float, line.split()[-3:]))
                forces_list.append(forces)

            if "Energy" in line:
                line = f.readline()
                energy = float(line.split()[0]) / natoms
                energy_list.append(energy)

            if "PlusStress" in line:
                line = f.readline()
                # MTP cfg stress 分量顺序
                plus_stress_mtp = np.array(list(map(float, line.split())))
                # xyz stress 分量顺序
                stress = [
                    plus_stress_mtp[0],
                    plus_stress_mtp[1],
                    plus_stress_mtp[2],
                    plus_stress_mtp[5],
                    plus_stress_mtp[3],
                    plus_stress_mtp[4],
                ]
                stress_list.append(stress)

                flag += 1
                print(f"No. {flag} configuration processed!")

    energy_array = np.array(energy_list).reshape(-1, 1)
    forces_array = np.concatenate(forces_list, axis=0)
    stress_array = np.array(stress_list)

    return energy_array, forces_array, stress_array


def main():
    dft_cfg_fn = "train.cfg"
    mtp_cfg_fn = "train_predict.cfg"

    dft_energy, dft_forces, dft_stress = extract_cfg(dft_cfg_fn)
    mtp_energy, mtp_forces, mtp_stress = extract_cfg(mtp_cfg_fn)

    energy_total = np.concatenate((mtp_energy, dft_energy), axis=1)
    forces_total = np.concatenate((mtp_forces, dft_forces), axis=1)
    stress_total = np.concatenate((mtp_stress, dft_stress), axis=1)

    energy_train_fn = "energy_train.out"
    force_train_fn = "force_train.out"
    stress_train_fn = "stress_train.out"
    np.savetxt(energy_train_fn, energy_total, fmt="%.10f")
    np.savetxt(force_train_fn, forces_total, fmt="%.10f")
    np.savetxt(stress_train_fn, stress_total, fmt="%.10f")

    print("\nEnergy, Forces, Stress data shape:")

    print(
        energy_total.shape,
        forces_total.shape,
        stress_total.shape,
    )

    print("\nWork is done!")


if __name__ == "__main__":
    main()
