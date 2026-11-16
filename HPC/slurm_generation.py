#!/usr/bin/env python3

"""生成不同平台的 VASP/LAMMPS/Python/Bash 任务 Slurm 提交脚本"""

from __future__ import annotations

import argparse
from typing import Literal

# 超算平台分区
HPC_PARTITION = {"sy": "64c512g", "pi": "cpu"}

# 课题组服务器：分区 -> 节点（见 sinfo），CPU 任务只用 master/node1/node2
GROUP_PARTITION_NODES = {
    "cpu": ["node1", "node3", "node4"],
    "cpu1": ["master"],
    "cpu2": ["node2"],
}
GROUP_NODES = ["master", "node1", "node2"]
NODE_PARTITION = {n: p for p, nodes in GROUP_PARTITION_NODES.items() for n in nodes}

INTEL_MODULES = {
    "vasp": ["intel-oneapi-compilers/2021.4.0", "intel-oneapi-mpi/2021.4.0", "intel-oneapi-mkl/2021.4.0"],
    "lammps": [
        "intel-oneapi-compilers/2021.4.0",
        "intel-oneapi-mpi/2021.12.1",
        "intel-oneapi-mkl/2021.4.0",
        "intel-oneapi-tbb/2021.4.0",
    ],
}


def group_node_directives(nodes: list[str]) -> list[str]:
    """课题组服务器：按所选节点生成 -p/-w/-x，分区需与节点对应"""
    partitions = sorted({NODE_PARTITION[n] for n in nodes})
    lines = [f"#SBATCH -p {','.join(partitions)}"]
    if len(nodes) == 1:
        lines.append(f"#SBATCH -w {nodes[0]}")
    else:
        # -w 多节点意味着全部分配，因此改用排除未选节点
        exclude = [n for p in partitions for n in GROUP_PARTITION_NODES[p] if n not in nodes]
        if exclude:
            lines.append(f"#SBATCH -x {','.join(exclude)}")
    lines.append("#SBATCH --no-requeue")
    return lines


def slurm_generation(
    platform: Literal["sy", "pi", "master"] = "sy",
    calculation_type: Literal["vasp", "lammps", "misc"] = "vasp",
    num_cpus: int = 1,
    time: str = "72:00:00",
    nodes: list[str] | None = None,
):
    """生成不同平台的 VASP/LAMMPS/Python/Bash 任务 Slurm 提交脚本"""

    on_group = platform == "master"
    # sy 平台未满 64 核时用 --ntasks 以便与他人共享节点
    ntasks = "ntasks" if platform == "sy" and num_cpus < 64 else "ntasks-per-node"

    lines = ["#!/bin/bash", "", f"#SBATCH -J {calculation_type.upper()}"]
    if on_group:
        lines += group_node_directives(nodes or ["node2"])
    else:
        lines.append(f"#SBATCH -p {HPC_PARTITION[platform]}")
    lines += [
        "#SBATCH -N 1",
        f"#SBATCH --{ntasks}={num_cpus}",
        f"#SBATCH -t {time}",
        "#SBATCH -o %j.out",
        "#SBATCH -e %j.err",
        "",
    ]

    if calculation_type in INTEL_MODULES and not on_group:
        modules = INTEL_MODULES[calculation_type]
        if calculation_type == "lammps" and platform == "sy":
            modules = modules + ["gcc/12.3.0"]
        lines += ["module purge"] + [f"module load {m}" for m in modules] + [""]

    bin_prefix = "" if on_group else "${HOME}/yangsl/bin/"
    if calculation_type == "vasp":
        lines += ["ulimit -s unlimited", "ulimit -l unlimited", "", f"mpirun {bin_prefix}vasp.5.std"]
    elif calculation_type == "lammps":
        lines.append(f"mpirun {bin_prefix}lmp_all -i in.lmp -sf opt -nocite")
    else:
        lines.append("# Please add your command here.")

    with open("job.slurm", "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"\n{platform} platform {calculation_type} submission file generated.")


def parse_nodes(value: str) -> list[str]:
    nodes = [n.strip() for n in value.split(",") if n.strip()]
    bad = [n for n in nodes if n not in GROUP_NODES]
    if not nodes or bad:
        raise argparse.ArgumentTypeError(f"invalid node(s) {bad or value!r}, choose from {GROUP_NODES}")
    return list(dict.fromkeys(nodes))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate slurm submission file for CPU jobs in group server and HPC platforms (VASP/LAMMPS/Misc).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument("platform", nargs="?", default="sy", choices=["sy", "pi", "master"], help="platform")
    parser.add_argument(
        "-ct", "--calculation_type", default="vasp", type=str.lower,
        choices=["vasp", "lammps", "misc"], metavar="STR", help="calculation type",
    )
    parser.add_argument("-nc", "--num_cpus", type=int, default=1, metavar="N", help="number of cpus")
    parser.add_argument("-t", "--time", default="72:00:00", metavar="STR", help="wall time")
    parser.add_argument(
        "--nodes", type=parse_nodes, default="node2", metavar="STR",
        help=f"node(s) for master platform from {GROUP_NODES}, e.g. node2 or master,node1,node2",
    )
    args = parser.parse_args()

    slurm_generation(
        platform=args.platform,
        calculation_type=args.calculation_type,
        num_cpus=args.num_cpus,
        time=args.time,
        nodes=args.nodes,
    )
