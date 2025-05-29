#!/usr/bin/env python3

"""生成 GPU 任务的 slurm 提交脚本（课题组服务器平台）"""

import argparse
from typing import Literal


# [ ] 添加 VASP GPU、LAMMPS GPU 任务
def write_slurm(
    num_cpus: int = 1,
    num_gpus: Literal[1, 2] = 1,
    calculation_type: Literal["nep", "gpumd"] = "gpumd",
):
    """生成 GPU 任务的 slurm 提交脚本（课题组服务器平台）"""

    with open("job.slurm", "w") as f:
        f.write(
            f"""#!/bin/bash

#SBATCH -J GPU
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks-per-node={num_cpus}
#SBATCH -t 72:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err

#SBATCH --gres=gpu:{num_gpus}

#SBATCH --no-requeue


if [[ "$SLURMD_NODENAME" == 'master' ]]; then
  {calculation_type}
elif [[ "$SLURMD_NODENAME" == 'node2' ]]; then
  {calculation_type}_node2
fi
"""
        )

    print(f"\n{calculation_type} submission file generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate slurm submission file for GPU jobs in group server platform (NEP/GPUMD).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "-nc",
        "--num_cpus",
        nargs="?",
        const=1,
        default=1,
        type=int,
        help="number of cpus, 1 or 2 is enough",
    )

    parser.add_argument(
        "-ng",
        "--num_gpus",
        nargs="?",
        const=1,
        default=1,
        type=int,
        choices=[1, 2],
        help="number of gpus",
    )

    parser.add_argument(
        "-ct",
        "--calculation_type",
        choices=["nep", "gpumd"],
        const="gpumd",
        default="gpumd",
        nargs="?",
        help="calculation type",
    )

    args = parser.parse_args()
    num_cpus = args.num_cpus
    num_gpus = args.num_gpus
    calculation_type = args.calculation_type

    write_slurm(num_cpus, num_gpus, calculation_type)
