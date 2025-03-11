"""
提取 OUTCAR 文件中的每个离子步中的每个原子受力数据并绘制

reference: https://github.com/shera-amit/utils/blob/main/scripts/bash_utils/plot_force
"""

import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from spt.plot_params import set_roman_plot_params

# [ ] 添加提取能量并绘制的功能


def extract_forces(outcar_fn: str = "OUTCAR") -> pd.DataFrame:
    """提取 OUTCAR 文件中的每个离子步中的每个原子受力数据"""

    with open(outcar_fn, "r") as f:
        lines = f.readlines()

    forces_pattern = re.compile(r"POSITION")
    forces_step_list = []
    forces_atom_list = []
    for i, line in enumerate(lines):
        if forces_pattern.search(line):
            j = i + 1
            while True:
                j += 1
                if (
                    lines[j].strip()
                    == "-----------------------------------------------------------------------------------"
                ):
                    forces_step_list.append(forces_atom_list)
                    forces_atom_list = []
                    break
                else:
                    force_vector = list(map(float, lines[j].strip().split()[-3:]))
                    forces_atom_list.append(np.linalg.norm(force_vector))

    forces_df = pd.DataFrame(forces_step_list).transpose()

    return forces_df


def plot_forces(forces_df: pd.DataFrame):
    """绘制最大、最小、平均力随离子步的变化"""

    min_force = forces_df.min(axis=0)
    max_force = forces_df.max(axis=0)
    average_force = forces_df.mean(axis=0)

    set_roman_plot_params()
    fig, ax = plt.subplots(figsize=(10, 6))

    x = list(range(1, forces_df.shape[1] + 1))
    ax.plot(x, average_force, marker="o", label="Average")
    ax.plot(x, min_force, marker="o", label="Minimum")
    ax.plot(x, max_force, marker="o", label="Maximum")

    ax.legend()

    ax.set_xlabel("Ion Step")
    ax.set_ylabel("Force (eV/Å)")

    ax.grid(True)

    fig.savefig("forces_step.png")


def main():
    outcar_fn = sys.argv[1]
    if not os.path.exists(outcar_fn):
        print("No OUTCAR file found.")
    else:
        forces_df = extract_forces(outcar_fn)
        plot_forces(forces_df)


if __name__ == "__main__":
    main()
