#!/usr/bin/env python3

"""二元 0K 计算相图绘制"""

import argparse
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator
from spt.plot_params import set_roman_plot_params


def to_subscript(formula_list: list[str]) -> list[str]:
    """将原化学式转换下标格式的字符串，供后续 LaTeX 渲染"""

    formula_subscript_list = [
        re.sub(r"([A-Z][a-z]*)(\d+)", r"\1_{\2}", formula) for formula in formula_list
    ]
    return formula_subscript_list


def binary_pd_plot(data_fn: str, element_list: list[str]):
    """二元 0K 计算相图绘制"""

    df = pd.read_csv(data_fn)

    matrix, solute = element_list

    # 按照 Ti-X 中的 X 元素含量值进行从小到大排序
    df.sort_values(by=[solute], inplace=True)
    # stable_df 用于绘制 Hull
    stable_df = df[df["e_above_hull"] == 0.0]
    # unstable_df 用于绘制散点图
    unstable_df = df[(df["e_above_hull"] > 0.0) & (df["e_above_hull"] < 0.1)]

    set_roman_plot_params()
    fig, ax = plt.subplots(figsize=(6, 6))

    # 散点数据绘制
    ax.scatter(unstable_df[solute], unstable_df["fepa"], s=50)

    stable_formula_list = stable_df["formula"]
    stable_formula_subscript_list = to_subscript(stable_formula_list)

    # marker 标记稳定相
    markers = ["p", "s", "v", "^", "D", "X", ">", "<"]
    for i in range(stable_df.shape[0]):
        ax.scatter(
            stable_df[solute].iloc[i],
            stable_df["fepa"].iloc[i],
            marker=markers[i],
            label=rf"$\rm{{{stable_formula_subscript_list[i]}}}$",
        )

    # Hull 绘制
    ax.plot(
        stable_df[solute],
        stable_df["fepa"],
        "-",
        color="green",
        # label="Hull",
    )

    x = np.arange(0, 1.1, 0.2)
    ax.set_xticks(x, labels=[f"{i:.1f}" for i in x])
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.set_ylim(-1.0, 0.2)
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.3))

    ax.set_xlabel(f"Composition ({solute})")
    ax.set_ylabel("Formation Energy (eV/atom)")

    fig_fn = f"{'_'.join(element_list)}.png"
    fig.savefig(fig_fn)

    print(f"Figure of {' '.join(element_list)} binary 0K phase diagram is generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot binary 0K phase diagram.",
        epilog="Author: SLY.",
    )

    parser.add_argument("data_fn", help="data file")
    parser.add_argument("element_list", nargs=2, help="element list (e.g. Ti Al)")

    args = parser.parse_args()

    binary_pd_plot(data_fn=args.data_fn, element_list=args.element_list)
