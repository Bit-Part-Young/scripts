#!/usr/bin/env python3

"""生成 GPU 任务的 Slurm 提交脚本（课题组服务器平台）"""

from __future__ import annotations

import argparse
from typing import Literal

# gpu 分区节点及其 GPU 数
GPU_NODES = {"master": 2, "node2": 1}


# [ ] 添加 VASP GPU、LAMMPS GPU 任务
def write_slurm(
    num_cpus: int = 1,
    num_gpus: int = 1,
    gpu_index: int | None = None,
    platform: Literal["master", "node2"] = "node2",
    calculation_type: Literal["nep", "gpumd"] = "gpumd",
    time: str = "72:00:00",
):
    """生成 GPU 任务的 Slurm 提交脚本（课题组服务器平台）"""

    if num_gpus > GPU_NODES[platform]:
        raise SystemExit(f"{platform} has only {GPU_NODES[platform]} GPU(s).")

    # 多卡时每块 GPU 至少配 1 个 CPU 核
    num_cpus = max(num_cpus, num_gpus)

    # 仅 master 单卡时可手动指定 GPU 编号
    cmd = calculation_type
    if platform == "master" and num_gpus == 1 and gpu_index is not None:
        cmd = f"CUDA_VISIBLE_DEVICES={gpu_index} {cmd}"

    with open("job.slurm", "w") as f:
        f.write(
            f"""#!/bin/bash

#SBATCH -J {calculation_type.upper()}
#SBATCH -p gpu
#SBATCH -w {platform}
#SBATCH -N 1
#SBATCH --ntasks-per-node={num_cpus}
#SBATCH --gres=gpu:{num_gpus}
#SBATCH -t {time}
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --no-requeue

{cmd}
"""
        )

    print(f"\n{platform} {calculation_type} submission file generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Slurm submission file for GPU jobs in group server platform (NEP/GPUMD).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument("-p", "--platform", choices=list(GPU_NODES), default="node2", metavar="STR",
                        help=f"node, one of {list(GPU_NODES)}")
    parser.add_argument("-nc", "--num_cpus", type=int, default=1, metavar="N", help="number of cpus, 1 or 2 is enough (at least num_gpus)")
    parser.add_argument("-ng", "--num_gpus", type=int, choices=[1, 2], default=1, metavar="N",
                        help="number of gpus (2 only on master)")
    parser.add_argument("-gi", "--gpu_index", type=int, choices=[0, 1], default=None, metavar="N",
                        help="gpu index on master (single gpu only); default lets Slurm choose")
    parser.add_argument("-ct", "--calculation_type", choices=["nep", "gpumd"], default="gpumd", metavar="STR",
                        help="calculation type")
    parser.add_argument("-t", "--time", default="72:00:00", metavar="STR", help="wall time")
    args = parser.parse_args()

    write_slurm(
        num_cpus=args.num_cpus,
        num_gpus=args.num_gpus,
        gpu_index=args.gpu_index,
        platform=args.platform,
        calculation_type=args.calculation_type,
        time=args.time,
    )
