#!/usr/bin/env python3

"""
NEP 势函数声子谱计算

Original Author: Ke Xu
"""

import argparse
import os

import numpy as np
import pandas as pd
from ase.atoms import Atoms
from ase.build import bulk
from ase.constraints import ExpCellFilter
from ase.io import read, write
from ase.optimize import BFGS, FIRE, GPMin
from ase.units import GPa
from calorine.tools import get_force_constants, relax_structure
from matplotlib import pyplot as plt
from seekpath import get_explicit_k_path

from phonopy.units import THzToCm


def Get_Disp(
    primitive_cell: Atoms,
    model_fn,
    potential_type="nep",
    nrep=2,
    fout="disp",
    Met_list=["C"],
    custom_path=None,
):
    """
    计算声子色散关系主函数

    params：
        unit_cell: ASE原子结构对象
        fpot: 势函数文件路径/赝势
        ptype: 势函数类型 nep/eam/tersoff/vasp
        nrep: 超胞扩展倍数
        fout: 输出文件名前缀
        Met_list: 元素列表
        custom_path: 自定义 k 点路径
    """

    # 势函数选择与结构优化
    if potential_type == "nep":
        from calorine.calculators import CPUNEP

        calculator = CPUNEP(model_fn)
        primitive_cell.calc = calculator
        relax_structure(
            primitive_cell,
            fmax=0.00001,
            minimizer="bfgs",
            constant_volume=False,
            constant_cell=False,
        )
    elif potential_type == "eam":
        from ase.calculators.lammpsrun import LAMMPS

        met = Met_list.join(" ")
        # print(met)
        parameters = {
            "pair_style": "eam/alloy",
            "pair_coeff": [f"* * {model_fn} {met}"],
        }
        calculator = LAMMPS(parameters=parameters, files=[model_fn])
        primitive_cell.set_calculator(calculator)
        mask = [True, True, True, True, True, True]
        ucf = ExpCellFilter(
            primitive_cell, scalar_pressure=0 * GPa, mask=mask, constant_volume=True
        )
        gopt = BFGS(ucf, maxstep=0.05)
        gopt.run(fmax=0.00001, steps=1000)
    elif potential_type == "tersoff":
        from ase.calculators.lammpsrun import LAMMPS

        mets = " ".join(Met_list)
        parameters = {
            "pair_style": "tersoff",
            "pair_coeff": ["* * {} {}".format(model_fn, mets)],
        }
        calculator = LAMMPS(files=[model_fn], **parameters)
        primitive_cell.set_calculator(calculator)
        mask = [True, True, True, True, True, True]
        ucf = ExpCellFilter(
            primitive_cell, scalar_pressure=0 * GPa, mask=mask, constant_volume=True
        )
        gopt = BFGS(ucf, maxstep=0.05)
        gopt.run(fmax=0.00001, steps=1000)
    elif potential_type == "vasp":
        from ase.calculators.vasp import Vasp

        calculator = Vasp(
            label="mylabel",
            txt="vasp.out",
            directory="VASP",
            command="/path/vasp_std",
            pp=model_fn,
        )
    else:
        raise "Only 'nep', 'eam', 'vasp' ptype could be used here."

    # k 点路径生成逻辑
    if custom_path is not None:
        # 自定义 k 点路径插值计算
        explicit_kpoints_rel = []
        explicit_labels = []
        num_points_per_segment = 68

        for start_label, end_label in custom_path:
            # 线性插值生成 k 点路径
            start = np.array(point_coords[start_label])
            end = np.array(point_coords[end_label])
            kpoints = [
                start + (end - start) * i / (num_points_per_segment - 1)
                for i in range(num_points_per_segment)
            ]
            explicit_kpoints_rel.extend(kpoints)
            # 路径标签处理
            if not explicit_labels or explicit_labels[-1] != start_label:
                explicit_labels.append(start_label)
            explicit_labels.extend([""] * (num_points_per_segment - 2))
            explicit_labels.append(end_label)
        explicit_kpoints_rel = np.array(explicit_kpoints_rel)
        # Create linear coordinates for plotting
        explicit_kpoints_linearcoord = np.linspace(0, 1, len(explicit_kpoints_rel))
        path = {
            "explicit_kpoints_rel": explicit_kpoints_rel,
            "explicit_kpoints_labels": explicit_labels,
            "explicit_kpoints_linearcoord": explicit_kpoints_linearcoord,
        }
    else:
        structure_tuple = (
            primitive_cell.cell,
            primitive_cell.get_scaled_positions(),
            primitive_cell.numbers,
        )
        path = get_explicit_k_path(structure_tuple, reference_distance=0.02)
        point_coords = path["point_coords"]

    # 声子计算与结果保存
    phonon = get_force_constants(primitive_cell, calculator, [nrep, nrep, nrep])
    phonon.run_band_structure([path["explicit_kpoints_rel"]])
    band = phonon.get_band_structure_dict()

    # 数据保存
    df = pd.DataFrame(band["frequencies"][0])
    np.savetxt(fout + ".dat", np.c_[df.index, df])
    np.save(fout + ".npy", path, allow_pickle=True)

    return df, path


def Get_Bandconf(
    path,
    atom_name="metals",
    nrep=3,
    ptype="nep",
    fout="band.conf",
):
    """
    生成Phonopy配置文件
    参数：
        path: k点路径数据
        atom_name: 元素名称
        nrep: 超胞扩展倍数
        ptype: 势函数类型
        fout: 输出文件名
    """

    # k点标签格式化处理
    kpoints = "BAND_LABELS = "
    coord = "BAND = "
    for i, pi in enumerate(path["path"]):
        # 特殊符号处理（Γ点）
        label = "$\Gamma$" if pi[0] == "GAMMA" else pi[0]
        kpoints += f"{label} "
        # 坐标转换
        coord += " ".join(map(str, path["point_coords"][pi[0]])) + " "

    # 配置文件写入
    with open(fout, "w") as fo:
        config = f"ATOM_NAME = {atom_name}\nDIM = {nrep} {nrep} {nrep}\n"
        config += f"{kpoints}\n{coord}\nFORCE_CONSTANTS = READ\n"
        fo.write(config)


def phonon_plot(df, path, pout="disp"):
    """
    声子色散曲线绘制函数
    参数：
        df: 频率数据
        path: k点路径
        pout: 输出图片前缀
    """

    # 绘图参数设置
    plt.figure(figsize=(4.2, 3), dpi=140)
    ax = plt.gca()

    # k点标签美化处理
    labels = [
        "$\Gamma$" if m == "GAMMA" else m.replace("_", "$_") + "$"
        for m in path["explicit_kpoints_labels"]
    ]

    # 频率曲线绘制
    for col in df.columns:
        ax.plot(df.index, df[col], color="cornflowerblue", lw=1)

    # 双坐标轴设置
    ax.set_ylabel("Frequency (THz)")
    ax2 = ax.twinx()
    ax2.set_ylabel("Frequency (cm$^{-1}$)")
    ax2.set_ylim(THzToCm * np.array(ax.get_ylim()))

    plt.tight_layout()
    plt.savefig(pout + ".png", dpi=300)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Phonon calculation using NEP potential."
    )
    parser.add_argument("structure_fn", type=str, help="structure filename")

    args = parser.parse_args()

    nrep = 4

    unit_cell = read(args.structure_fn)
    cal_type = "nep"
    met = "".join(list(set(unit_cell.get_chemical_symbols())))

    # 计算目录管理
    folder = f"{met}-{cal_type}"
    os.makedirs(folder, exist_ok=True)

    # 执行计算流程
    df, path = Get_Disp(
        unit_cell,
        model_fn="nep.txt",
        potential_type=cal_type,
        nrep=nrep,
        fout=f"{folder}/disp.dat",
    )
    Get_Bandconf(path, atom_name=met, nrep=nrep, fout=f"{folder}/band.conf")
    phonon_plot(df, path, pout=f"{folder}/disp_plot")
