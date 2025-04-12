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
    calculation_type: Literal["VASP", "LAMMPS", "Python", "Bash"] = "VASP",
    slurm_fn: str = "job.slurm",
    num_cpus: int = 1,
    time: str = "72:00:00",
    input_fn: str = "in.lmp",
):
    """生成不同平台的 VASP/LAMMPS/Python/Bash 任务 Slurm 提交脚本"""

    with open(slurm_fn, "w") as f:
        f.write("#!/bin/bash\n\n")

        partition = platform_partition[platform]

        f.write(f"#SBATCH -J {calculation_type}\n")
        f.write(f"#SBATCH -p {partition}\n")
        f.write("#SBATCH -N 1\n")
        f.write(f"#SBATCH --ntasks-per-node={num_cpus}\n")
        f.write(f"#SBATCH -t {time}\n")

        f.write("#SBATCH -o %j.out\n")
        f.write("#SBATCH -e %j.err\n\n")

        if platform == "master":
            f.write("#SBATCH -w node2\n")
            f.write("#SBATCH -x node1\n\n")
            f.write("#SBATCH --no-requeue\n\n")

        if calculation_type == "VASP":
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

        elif calculation_type == "LAMMPS":
            if platform == "sy":
                lammps_cmd = "${HOME}/yangsl/bin/lmp_all"
            elif platform == "master":
                lammps_cmd = "lmp_all"

                f.write("module load intel-oneapi-compilers/2021.4.0\n")
                f.write("module load intel-oneapi-mpi/2021.12.1\n")
                f.write("module load intel-oneapi-mkl/2021.4.0\n")
                f.write("module load intel-oneapi-tbb/2021.4.0\n\n")

            f.write(f"mpirun {lammps_cmd} -i {input_fn}\n")

        elif calculation_type == "Python":
            f.write(f"conda activate base_ysl\n\n")

            f.write(f"python {input_fn}\n")

        elif calculation_type == "Bash":
            f.write(f"bash {input_fn}\n")

    print(f"\n{platform} platform {calculation_type} submission file generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate submission file for VASP/LAMMPS/Python/Bash calculation in HPC.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "platform",
        nargs="?",
        const="sy",
        default="sy",
        type=str,
        choices=["sy", "pi", "master"],
        help="Platform",
    )

    parser.add_argument(
        "-ct",
        "--calculation_type",
        nargs="?",
        const="VASP",
        default="VASP",
        type=str,
        choices=["VASP", "LAMMPS", "Python", "Bash"],
        help="Calculation type",
    )

    parser.add_argument(
        "-nc",
        "--num_cpus",
        nargs="?",
        const=1,
        default=1,
        type=int,
        help="Number of CPUs",
    )

    parser.add_argument(
        "--input_fn",
        nargs="?",
        const="in.lmp",
        default="in.lmp",
        type=str,
        help="LAMMPS/Python/Bash input filename",
    )

    args = parser.parse_args()

    slurm_generation(
        platform=args.platform,
        calculation_type=args.calculation_type,
        num_cpus=args.num_cpus,
        input_fn=args.input_fn,
    )
