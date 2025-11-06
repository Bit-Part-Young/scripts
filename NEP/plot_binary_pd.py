#!/usr/bin/env python3

"""使用 pymatgen 绘制 0K 二元相图（根据数据文件中的 formula 和 energy 数据列绘制）"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry

from spt.plot_params import set_plot_params


def plot_binary_pd(data_fn: str, cols: list[int], fig_fn: str):
    """使用 pymatgen 绘制 0K 二元相图（根据数据文件中的 formula 和 energy 数据列绘制）"""

    df = pd.read_csv(data_fn, sep=None, engine="python")

    entries: list[PDEntry] = df.apply(
        lambda row: PDEntry(
            composition=Composition(row.iloc[cols[0]]), energy=row.iloc[cols[1]]
        ),
        axis=1,
    ).tolist()

    phase_diagram = PhaseDiagram(entries)
    # stable_entries 类型为 set
    stable_entries = phase_diagram.stable_entries
    print(f"\nStable entries info:")
    for stable_entry in stable_entries:
        print(stable_entry.composition, stable_entry.energy)
    print("\nNote: The energy in MP is not real DFT calculated energy?")

    set_plot_params(roman_params=True, sci_params=True)

    # [ ] 绘图参数待优化
    ax = phase_diagram.get_plot(
        backend="matplotlib",  # 绘图后端；ploty 或 matplotlib
        show_unstable=True,  # 是否显示非稳定构型
        # label_stable=False,  # 是否显示 label
    )

    plt.savefig(fig_fn)

    print(f"\n0K binary phase diagram {fig_fn} is generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot binary 0K phase diagram with pymatgen from formula and energy data",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "data_fn", help="data filename (must contain formula and energy columns)"
    )
    parser.add_argument(
        "cols",
        type=int,
        nargs=2,
        help="column numbers for formula and energy (starting from 0)",
    )
    parser.add_argument(
        "-o",
        "--fig_fn",
        default="binary_pd.png",
        metavar="FILE",
        help="output figure filename",
    )

    args = parser.parse_args()

    plot_binary_pd(args.data_fn, args.cols, args.fig_fn)
