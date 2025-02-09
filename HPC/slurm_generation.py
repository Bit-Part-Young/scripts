"""生成 VASP/LAMMPS/Python/Bash 任务 slurm 提交脚本"""

import argparse


def slurm_generation(
    calculation_type: str = "VASP",
    slurm_file: str = "job.slurm",
    partition: str = "64c512g",
    num_cpus: int = 1,
    time: str = "72:00:00",
    input_fn: str = "in.lmp",
):
    """生成 VASP/LAMMPS/Python/Bash 任务 slurm 提交脚本"""

    with open(slurm_file, "w") as f:
        f.write("#!/bin/bash\n\n")

        f.write(f"#SBATCH -J {calculation_type}\n")
        f.write(f"#SBATCH -p {partition}\n")
        f.write("#SBATCH -N 1\n")
        f.write(f"#SBATCH --ntasks-per-task={num_cpus}\n")
        f.write(f"#SBATCH -t {time}\n")
        f.write("#SBATCH -o %j.out\n")
        f.write("#SBATCH -e %j.err\n\n")

        f.write("module purge\n\n")

        if calculation_type == "VASP":
            vasp_cmd = "${HOME}/yangsl/bin/vasp.5.std"

            f.write("module load intel-oneapi-compilers/2021.4.0\n")
            f.write("module load intel-oneapi-mpi/2021.4.0\n")
            f.write("module load intel-oneapi-mkl/2021.4.0\n\n")

            f.write("ulimit -s unlimited\n")
            f.write("ulimit -l unlimited\n\n")

            f.write(f"mpirun {vasp_cmd}\n")

        elif calculation_type == "LAMMPS":

            lammps_cmd = "${HOME}/yangsl/bin/lmp_all"

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

    print(f"\n{calculation_type} submission file generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate submission file for VASP/LAMMPS/Python/Bash calculation in HPC.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        "-p",
        "--partition",
        nargs="?",
        const="64c512g",
        default="64c512g",
        choices=["64c512g", "cpu"],
        type=str,
        help="Partition name",
    )

    parser.add_argument(
        "--input_fn",
        nargs="?",
        const="in.lmp",
        default="in.lmp",
        type=str,
        help="Input filename",
    )

    args = parser.parse_args()

    slurm_generation(
        calculation_type=args.calculation_type,
        num_cpus=args.num_cpus,
        partition=args.partition,
        input_fn=args.input_fn,
    )
