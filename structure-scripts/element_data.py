#!/usr/bin/env python3

"""
元素周期表常见元素数据

atomic_mass: 相对原子质量
atomic_radius: 原子半径（无统一定值，不一定准确）
melting_point: 熔点（K）
density: 密度（g/cm^3）
energy_pa: 基态结构的平均原子能量（eV/atom）
"""

import argparse

import pandas as pd

data = {
    "Ti": (47.867, 1.4, 1941.0, 4.54, -7.8),
    "Al": (26.982, 1.25, 933.47, 2.7, -3.7),
    "Nb": (92.906, 1.45, 2750.0, 8.57, -10.2),
    "Zr": (91.224, 1.55, 2128.0, 6.52, -8.5),
    "V": (50.941, 1.35, 2183.0, 6.11, -8.9),
    "Mo": (95.94, 1.45, 2896.0, 10.28, -10.9),
    "Ta": (180.948, 1.45, 3290.0, None, -11.8),
    "W": (183.84, 1.35, 3695.0, None, -13.0),
    "Hf": (178.49, 1.55, 2506.0, None, -9.9),
    "Sn": (118.71, 1.45, 505.08),
    "Mg": (24.305, 1.5, 923.0),
    "Cu": (63.546, 1.35, 1357.77),
    "Ni": (58.693, 1.35, 1728.0),
    "Cr": (51.996, 1.4, 2180.0),
    "Ag": (107.868, 1.6, 1234.93),
    "Au": (196.967, 1.35, 1337.33),
    "Pb": (207.2, 1.8, 600.61),
    "Pd": (106.42, 1.4, 1828.05),
    "Pt": (195.084, 1.35, 2041.4),
}

df = pd.DataFrame.from_dict(
    data,
    orient="index",
    columns=[
        "atomic_mass",
        "atomic_radius",
        "melting_point(K)",
        "density(g/cm^3)",
        "energy_pa(eV/atom)",
    ],
).fillna("-")


parser = argparse.ArgumentParser(
    description="Common elements data.",
    # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    epilog="Author: SLY.",
)

parser.add_argument(
    "-e",
    "--element",
    type=str,
    nargs="+",
    metavar="element",
    help="Element symbol, e.g. Ti",
)

args = parser.parse_args()

if args.element:
    print(df.loc[args.element])
else:
    print(df)
