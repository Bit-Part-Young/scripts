#!/usr/bin/env python3

"""生成不同平台的 VASP/LAMMPS/Python/Bash 任务 Slurm 提交脚本"""

import argparse
from typing import Literal

platform_partition = {
    "sy": "64c512g",
    "pi": "cpu",
    "master": "cpu",
}


def slurm_generation(
    platform: Literal["sy", "pi", "master"] = "sy",
    calculation_type: Literal["vasp", "lammps", "misc"] = "vasp",
    num_cpus: int = 1,
    time: str = "72:00:00",
):
    """生成不同平台的 VASP/LAMMPS/Python/Bash 任务 Slurm 提交脚本"""

    with open("job.slurm", "w") as f:
        f.write("#!/bin/bash\n\n")

        f.write(f"#SBATCH -J {calculation_type.upper()}\n")

        partition = platform_partition[platform]
        f.write(f"#SBATCH -p {partition}\n")

        f.write("#SBATCH -N 1\n")
        f.write(f"#SBATCH --ntasks-per-node={num_cpus}\n")
        f.write(f"#SBATCH -t {time}\n")

        f.write("#SBATCH -o %j.out\n")
        f.write("#SBATCH -e %j.err\n\n")

        if platform == "master":
            f.write("#SBATCH -x master\n\n")
            f.write("#SBATCH --no-requeue\n\n")

        if calculation_type == "vasp":
            if platform == "master":
                vasp_cmd = "vasp.5.std"
            else:
                vasp_cmd = "${HOME}/yangsl/bin/vasp.5.std"

                f.write("module purge\n\n")
                f.write("module load intel-oneapi-compilers/2021.4.0\n")
                f.write("module load intel-oneapi-mpi/2021.4.0\n")
                f.write("module load intel-oneapi-mkl/2021.4.0\n\n")

            f.write("ulimit -s unlimited\n")
            f.write("ulimit -l unlimited\n\n")

            f.write(f"mpirun {vasp_cmd}\n")

        elif calculation_type == "lammps":
            if platform == "master":
                lammps_cmd = "lmp_all"
            else:
                lammps_cmd = "${HOME}/yangsl/bin/lmp_all"

                f.write("module load intel-oneapi-compilers/2021.4.0\n")
                f.write("module load intel-oneapi-mpi/2021.12.1\n")
                f.write("module load intel-oneapi-mkl/2021.4.0\n")
                f.write("module load intel-oneapi-tbb/2021.4.0\n\n")

            f.write(f"mpirun {lammps_cmd} -i in.lmp\n")

        elif calculation_type == "misc":
            f.write(f"# Please add your command here.\n\n")

    print(f"\n{platform} platform {calculation_type} submission file generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate slurm submission file for CPU jobs in group server and HPC platforms (VASP/LAMMPS/Misc).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "platform",
        nargs="?",
        const="sy",
        default="sy",
        choices=["sy", "pi", "master"],
        help="Platform",
    )

    parser.add_argument(
        "-ct",
        "--calculation_type",
        nargs="?",
        const="vasp",
        default="vasp",
        choices=["vasp", "lammps", "misc"],
        help="calculation type",
    )

    parser.add_argument(
        "-nc",
        "--num_cpus",
        nargs="?",
        const=1,
        default=1,
        type=int,
        help="number of cpus",
    )

    args = parser.parse_args()

    platform = args.platform
    calculation_type: str = args.calculation_type
    num_cpus = args.num_cpus

    slurm_generation(
        platform=platform,
        calculation_type=calculation_type.lower(),
        num_cpus=num_cpus,
    )
