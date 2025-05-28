#!/usr/bin/env python3

"""
绘制 GPUMD 输出的 thermo.out 文件数据

reference: https://github.com/zhyan0603/GPUMDkit/blob/main/Scripts/plt_scripts/plt_nep_thermo.py
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.ticker import ScalarFormatter

from spt.plot_params import set_plot_params


def calculate_angle(x, y):
    dot_product = np.einsum("ij,ij->i", x, y)
    norm_x = np.linalg.norm(x, axis=1)
    norm_y = np.linalg.norm(y, axis=1)
    angle_radians = np.arccos(dot_product / (norm_x * norm_y))
    return np.degrees(angle_radians)


def calculate_volume(a, b, c):
    volume = np.einsum("ij,ij->i", a, np.cross(b, c))
    return np.abs(volume)


# Determine dump_interval from run.in file
def get_dump_interval():
    timestep = 1.0  # Default
    dump_interval = 10  # Default

    if os.path.exists("run.in"):
        with open("run.in", "r") as file:
            for line in file:
                # Read timestep value
                if "time_step" in line:
                    try:
                        timestep = float(line.split()[1])  # timestep in fs
                    except (IndexError, ValueError):
                        pass
                # Read dump_thermo interval
                elif "dump_thermo" in line:
                    try:
                        dump_interval = int(
                            line.split()[1]
                        )  # number of timesteps between dumps
                        break
                    except (IndexError, ValueError):
                        pass

    # Calculate total time interval per dump in ps
    total_interval_ps = timestep * dump_interval / 1000.0

    return total_interval_ps


data = np.loadtxt("./thermo.out")

dump_interval_ps = get_dump_interval()
time = np.arange(0, len(data) * dump_interval_ps, dump_interval_ps)

# read data
temperature = data[:, 0]
kinetic_energy = data[:, 1]
potential_energy = data[:, 2]
pressure_x = data[:, 3]
pressure_y = data[:, 4]
pressure_z = data[:, 5]

ncols = data.shape[1]
if ncols == 12:
    box_length_x = data[:, 9]
    box_length_y = data[:, 10]
    box_length_z = data[:, 11]
    volume = box_length_x * box_length_y * box_length_z
elif ncols == 18:
    ax, ay, az = data[:, 9], data[:, 10], data[:, 11]
    bx, by, bz = data[:, 12], data[:, 13], data[:, 14]
    cx, cy, cz = data[:, 15], data[:, 16], data[:, 17]

    a_vectors = np.column_stack((ax, ay, az))
    b_vectors = np.column_stack((bx, by, bz))
    c_vectors = np.column_stack((cx, cy, cz))

    box_length_x = np.sqrt(ax**2 + ay**2 + az**2)
    box_length_y = np.sqrt(bx**2 + by**2 + bz**2)
    box_length_z = np.sqrt(cx**2 + cy**2 + cz**2)

    box_angle_alpha = calculate_angle(b_vectors, c_vectors)
    box_angle_beta = calculate_angle(c_vectors, a_vectors)
    box_angle_gamma = calculate_angle(a_vectors, b_vectors)

    volume = calculate_volume(a_vectors, b_vectors, c_vectors)
else:
    raise ValueError("Unsupported number of columns in thermo.out. Expected 12 or 18.")

# Calculate averages after 50% of simulation time
start_index = int(len(time) * 0.5)  # You can change it based on your need
avg_temperature = np.mean(temperature[start_index:])
avg_pressure_x = np.mean(pressure_x[start_index:])
avg_pressure_y = np.mean(pressure_y[start_index:])
avg_pressure_z = np.mean(pressure_z[start_index:])
avg_length_x = np.mean(box_length_x[start_index:])
avg_length_y = np.mean(box_length_y[start_index:])
avg_length_z = np.mean(box_length_z[start_index:])
avg_volume = np.mean(volume[start_index:]) / 1000  # Convert to x10^3 Å^3
if ncols == 18:
    avg_angle_alpha = np.mean(box_angle_alpha[start_index:])
    avg_angle_beta = np.mean(box_angle_beta[start_index:])
    avg_angle_gamma = np.mean(box_angle_gamma[start_index:])
    avg_ax = np.mean(ax[start_index:])
    avg_ay = np.mean(ay[start_index:])
    avg_az = np.mean(az[start_index:])
    avg_bx = np.mean(bx[start_index:])
    avg_by = np.mean(by[start_index:])
    avg_bz = np.mean(bz[start_index:])
    avg_cx = np.mean(cx[start_index:])
    avg_cy = np.mean(cy[start_index:])
    avg_cz = np.mean(cz[start_index:])

# Print average values
average_results = [
    f"+------------------------------------------+",
    f"| Average values after 50% simulation time |",
    f"| You can change it by the following line: |",
    f"|    start_index = int(len(time) * 0.5)    |",
    f"+------------------------------------------+",
    f"Temperature: {avg_temperature:.3f} K",
    f"Pressure X: {avg_pressure_x:.3f} GPa",
    f"Pressure Y: {avg_pressure_y:.3f} GPa",
    f"Pressure Z: {avg_pressure_z:.3f} GPa",
    f"Lattice Length X: {avg_length_x:.3f} Å",
    f"Lattice Length Y: {avg_length_y:.3f} Å",
    f"Lattice Length Z: {avg_length_z:.3f} Å",
    f"Volume: {avg_volume*1000:.3f} Å^3",
]
if ncols == 18:
    average_results.extend(
        [
            f"Angle Alpha: {avg_angle_alpha:.2f}°",
            f"Angle Beta: {avg_angle_beta:.2f}°",
            f"Angle Gamma: {avg_angle_gamma:.2f}°",
        ]
    )
    average_results.append(
        f"Average lattice matrix: {avg_ax:.3f}, {avg_ay:.3f}, {avg_az:.3f}, {avg_bx:.3f}, {avg_by:.3f}, {avg_bz:.3f}, {avg_cx:.3f}, {avg_cy:.3f}, {avg_cz:.3f}"
    )

print("\n".join(average_results))

# Save average values to a text file
with open("./average_results.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(average_results))


def plot_gpumd_thermo():

    set_plot_params(savefig_dpi=300)

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # 温度
    ax: Axes = axs[0, 0]
    ax.plot(time, temperature)
    ax.set_title("Temperature")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Temperature (K)")

    # 压力
    ax: Axes = axs[0, 1]
    ax.plot(time, pressure_x, label="Px")
    ax.plot(time, pressure_y, label="Py")
    ax.plot(time, pressure_z, label="Pz")
    ax.set_title("Pressure")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Pressure (GPa)")
    ax.legend()

    # 势能和动能
    # Determine y-axis limits dynamically
    pe_min, pe_max = np.min(potential_energy), np.max(potential_energy)
    pe_range = pe_max - pe_min
    pe_ylim_lower = pe_min - 0.6 * pe_range  # Double the range downward
    pe_ylim_upper = pe_max + 0.05 * pe_range  # Double the range upward
    ax: Axes = axs[0, 2]
    ax.set_title(r"$P_E$ vs $K_E$")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"Potential Energy (eV)", color="tab:orange")
    ax.plot(time, potential_energy, color="tab:orange")
    ax.set_ylim(pe_ylim_lower, pe_ylim_upper)  # Set extended range for PE
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.tick_params(axis="y", labelcolor="tab:orange")

    ke_min, ke_max = np.min(kinetic_energy), np.max(kinetic_energy)
    ke_range = ke_max - ke_min
    ke_ylim_lower = ke_min - 0.05 * ke_range  # Double the range downward
    ke_ylim_upper = ke_max + 0.6 * ke_range  # Double the range upward
    ax_kinetic: Axes = axs[0, 2].twinx()
    ax_kinetic.set_ylabel("Kinetic Energy (eV)", color="tab:green")
    ax_kinetic.plot(time, kinetic_energy, color="tab:green")
    ax_kinetic.set_ylim(ke_ylim_lower, ke_ylim_upper)  # Set extended range for KE
    ax_kinetic.tick_params(axis="y", labelcolor="tab:green")
    ax_kinetic.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax_kinetic.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # 晶格参数
    ax: Axes = axs[1, 0]
    ax.plot(time, box_length_x, label="Lx")
    ax.plot(time, box_length_y, label="Ly")
    ax.plot(time, box_length_z, label="Lz")
    ax.set_title("Lattice Parameters")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"Lattice Parameters ($\AA$)")
    ax.legend()

    # 体积
    ax: Axes = axs[1, 1]
    ax.plot(time, volume, label="Volume", color="tab:purple")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_title("Volume")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(r"Volume ($\AA^3$)")
    ax.legend()

    # 夹角 (仅适用于三斜晶系)
    if ncols == 18:
        ax: Axes = axs[1, 2]
        ax.plot(time, box_angle_alpha, label=r"$\alpha$")
        ax.plot(time, box_angle_beta, label=r"$\beta$")
        ax.plot(time, box_angle_gamma, label=r"$\gamma$")
        ax.set_title("Interaxial Angles")
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(r"Interaxial Angles ($\degree$)")
        ax.legend()

    plt.tight_layout()

    fig.savefig("thermo.png")


if __name__ == "__main__":
    plot_gpumd_thermo()
