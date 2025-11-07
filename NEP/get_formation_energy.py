#!/usr/bin/env python3

"""获取给定构型的形成能（用于后续的 0K 二元相图绘制）"""

import argparse
import os

import pandas as pd
from ase.formula import Formula
from ase.io import read

# 优先以元素 Al 作为 0K 二元相图的 x 轴变量
element_sequence_list = ["Al", "Ti", "Nb", "Mo", "Zr", "V"]
# element_sequence_list = ["Ti", "Al", "Nb", "Mo", "Zr", "V"]


def get_element_energy_dict(element_energy_fn: str | None = None):
    if element_energy_fn is None:
        # DFT 元素基态结构的平均原子能量
        element_energy_dict = {
            "Ti": -7.83539,
            "Al": -3.74433,
            "Nb": -10.21808,
            "Mo": -10.93079,
            "Zr": -8.52201,
            "V": -8.99063,
        }
    else:
        # 从文件中读取元素基态结构的平均原子能量
        df = pd.read_csv(element_energy_fn, sep=None, engine="python")
        element_energy_dict = df.set_index("element")["energy"].to_dict()

    return element_energy_dict


def get_energy(
    structure_fn: str,
    add_element_data: bool = False,
    element_energy_fn: str | None = None,
):
    """获取给定构型的化学式和能量"""

    element_energy_dict = get_element_energy_dict(element_energy_fn)

    input_fn_dirname = os.path.dirname(structure_fn)
    input_fn_basename = os.path.basename(structure_fn)
    input_format = input_fn_basename.split(".")[-1]

    formula_energy_list = []
    if input_format in ["vasp", "POSCAR"]:
        outcar_fn = os.path.join(input_fn_dirname, "OUTCAR")
        if os.path.exists(outcar_fn):
            atoms = read(outcar_fn, format="vasp-out")

            formula = atoms.get_chemical_formula()
            energy = atoms.get_potential_energy()

            formula_energy_list.append((formula, energy))
        else:
            raise FileNotFoundError("Error: OUTCAR not found!")

    elif input_format in ["xyz", "extxyz"]:
        atoms_list = read(structure_fn, index=":", format="extxyz")

        for atoms in atoms_list:
            formula = atoms.get_chemical_formula()
            energy = atoms.get_potential_energy()

            formula_energy_list.append((formula, energy))

        # 是否添加元素对应的数据
        if add_element_data:
            composition = Formula(formula_energy_list[0][0]).count()
            element_list = list(composition.keys())

            # 实现 Ti1Al0 Ti1Al0Nb0 的化学式样式，便于后续的 0K 二元相图绘制
            for element in element_list:
                formula_partial = []
                for e in element_list:
                    count = 1 if e == element else 0
                    formula_partial.append(f"{e}{count}")
                formula = "".join(formula_partial)
                formula_energy_list.append((formula, element_energy_dict[element]))

    return formula_energy_list


def get_formation_energy(
    structure_fn: str,
    add_element_data: bool = False,
    element_energy_fn: str | None = None,
):
    """获取给定构型的化学式和形成能"""

    element_energy_dict = get_element_energy_dict(element_energy_fn)
    formula_energy_list = get_energy(structure_fn, add_element_data, element_energy_fn)

    formula_formation_energy_list = []
    data1_list = []
    data2_list = []
    for formula, energy in formula_energy_list:
        composition = Formula(formula).count()

        natoms = sum(composition.values())
        element_list = list(composition.keys())
        # 按照给定的元素顺序排序
        element_list.sort(key=lambda x: element_sequence_list.index(x))

        composition_fractional = {
            element: composition.get(element, 0) / natoms for element in element_list
        }

        # 形成能计算
        formation_energy = energy
        for element, count in composition.items():
            formation_energy -= element_energy_dict[element] * count
        formation_energy_pa = round(formation_energy / natoms, 5)

        formula_formation_energy_list.append((formula, formation_energy_pa))

        data1_dict = {
            "formula": formula,
            "energy": energy,
            "fe": formation_energy_pa,
        }
        data1_list.append(data1_dict)

        data2_dict = composition_fractional
        data2_dict.update({"fe": formation_energy_pa})
        data2_list.append(data2_dict)

    # 分别保存成 化学式-能量-形成能 和 成分-形成能 两个数据文件
    # 化学式-能量-形成能 是使用 plot_binary_pd.py 脚本绘制 0K 二元相图需要的数据
    # 成分-形成能 是使用 pycxl 和 plot_binary_pd.py 脚本绘制 0K 二元相图需要的数据
    df1 = pd.DataFrame(data1_list).round(5)
    df2 = pd.DataFrame(data2_list).round(5)
    df1.to_csv("formula_formation_energy.dat", sep=" ", index=False)
    df2.to_csv("composition_formation_energy.dat", sep=" ", index=False)

    print("\n2 data files generated, please look and check!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get formation energy of given configuration(s), useful for 0K binary phase diagram plotting.",
        epilog="Author: SLY.",
    )

    parser.add_argument("structure_fn", help="structure filename; eg. POSCAR, *.xyz")
    parser.add_argument(
        "--add_element_data", action="store_true", help="whether to add element data"
    )
    parser.add_argument(
        "--element_energy_fn",
        metavar="FILE",
        help="element energy filename",
    )

    args = parser.parse_args()

    get_formation_energy(
        args.structure_fn, args.add_element_data, args.element_energy_fn
    )
