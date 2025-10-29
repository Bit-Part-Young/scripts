#!/usr/bin/env python3

"""生成 GPU 任务的 Slurm 提交脚本（课题组服务器平台）"""

import argparse
from typing import Literal


# [ ] 添加 VASP GPU、LAMMPS GPU 任务
def write_slurm(
    num_cpus: int = 1,
    num_gpus: int = 1,
    gpu_index: int = 1,
    platform: Literal["master", "node2"] = "node2",
    calculation_type: Literal["nep", "gpumd"] = "gpumd",
):
    """生成 GPU 任务的 Slurm 提交脚本（课题组服务器平台）"""

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

#SBATCH -w {platform}

#SBATCH --gres=gpu:{num_gpus}

#SBATCH --no-requeue


if [[ "$SLURMD_NODENAME" == 'master' ]]; then
  if [[ $SLURM_GPUS_ON_NODE == 1 ]]; then
    CUDA_VISIBLE_DEVICES={gpu_index} {calculation_type}
  else
    {calculation_type}
  fi
elif [[ "$SLURMD_NODENAME" == 'node2' ]]; then
  {calculation_type}_node2
fi
"""
        )

    print(f"\n{calculation_type} submission file generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Slurm submission file for GPU jobs in group server platform (NEP/GPUMD).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "-p",
        "--platform",
        choices=["master", "node2"],
        default="node2",
        metavar="STR",
        help="platform",
    )

    parser.add_argument(
        "-nc",
        "--num_cpus",
        type=int,
        default=1,
        metavar="N",
        help="number of cpus, 1 or 2 is enough",
    )

    parser.add_argument(
        "-ng",
        "--num_gpus",
        type=int,
        choices=[1, 2],
        default=1,
        metavar="N",
        help="number of gpus",
    )

    parser.add_argument(
        "-gi",
        "--gpu_index",
        type=int,
        default=1,
        metavar="N",
        help="gpu index",
    )

    parser.add_argument(
        "-ct",
        "--calculation_type",
        choices=["nep", "gpumd"],
        default="gpumd",
        metavar="STR",
        help="calculation type",
    )

    args = parser.parse_args()

    write_slurm(
        num_cpus=args.num_cpus,
        num_gpus=args.num_gpus,
        gpu_index=args.gpu_index,
        platform=args.platform,
        calculation_type=args.calculation_type,
    )
