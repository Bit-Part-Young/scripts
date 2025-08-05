#!/usr/bin/env python3

"""
元素周期表常见元素数据

atomic_number: 原子序数
atomic_mass: 相对原子质量
atomic_radius: 原子半径（无统一定值，不一定准确）
melting_point (K): 熔点
density (g/cm^3): 密度
a(Å): 晶格常数 a
c(Å): 晶格常数 c
energy_pa (eV/atom): 基态结构的平均原子能量
mp-id: Materials Project ID
"""

import argparse

import pandas as pd

data = {
    "Ti": (22, 47.867, 1.4, 1941.0, 4.51, -7.8, 2.93, 4.65, "mp-46"),
    "Al": (13, 26.982, 1.25, 933.47, 2.70, -3.7, 4.04, 4.04, "mp-134"),
    "Nb": (41, 92.906, 1.45, 2750.0, 8.57, -10.2, 3.31, 3.31, "mp-75"),
    "Zr": (40, 91.224, 1.55, 2128.0, 6.51, -8.5, 3.23, 5.15, "mp-131"),
    "V": (23, 50.941, 1.35, 2183.0, 6.11, -8.9, 2.98, 2.98, "mp-146"),
    "Mo": (42, 95.94, 1.45, 2896.0, 10.28, -10.9, 3.16, 3.16, "mp-129"),
    "Ta": (73, 180.948, 1.45, 3290.0, 16.65, -11.8, 3.32, 3.32, "mp-50"),
    "W": (74, 183.84, 1.35, 3695.0, 19.25, -13.0, 3.18, 3.18, "mp-91"),
    "Hf": (72, 178.49, 1.55, 2506.0, 13.31, -9.9, 3.20, 5.06, "mp-103"),
    "Sn": (50, 118.71, 1.45, 505.08, 7.31, -3.8, 6.65, 6.65, "mp-117"),
    "Mg": (12, 24.305, 1.5, 923.0, 1.74, -1.5, 3.20, 5.16, "mp-153"),
    "Cu": (29, 63.546, 1.35, 1357.77, 8.92, -3.7, 3.63, 3.63, "mp-30"),
    "Si": (14, 28.0855, 1.1, 1687.0, 2.33, -5.4, 5.47, 5.47, "mp-149"),
    "Ni": (28, 58.693, 1.35, 1728.0, 8.91, -5.5, 3.52, 3.52, "mp-23"),
    "Cr": (24, 51.996, 1.4, 2180.0, 7.14, -9.5, 2.84, 2.84, "mp-90"),
    "Ag": (47, 107.868, 1.6, 1234.93, 10.49, -2.7, 4.14, 4.14, "mp-124"),
    "Au": (79, 196.967, 1.35, 1337.33, 19.30, -3.2, 4.16, 4.16, "mp-81"),
    "Pb": (82, 207.2, 1.8, 600.61, 11.34, None, None, None, "mp-20483"),
    "Pd": (46, 106.42, 1.4, 1828.05, 12.02, -5.2, 3.94, 3.94, "mp-2"),
    "Pt": (78, 195.084, 1.35, 2041.4, 21.09, -6.1, 3.97, 3.97, "mp-126"),
}

df = pd.DataFrame.from_dict(
    data,
    orient="index",
    columns=[
        "atomic_number",
        "atomic_mass",
        "atomic_radius",
        "T_m(K)",
        "ρ(g/cm^3)",
        "E_pa(eV/atom)",
        "a(Å)",
        "c(Å)",
        "mp-id",
    ],
).fillna("-")


# 改成英文
def note():
    print("\nNote:")
    print(
        "1. atomic_radius: atomic radius (no unified value, not necessarily accurate), got from pymatgen"
    )
    print(
        "2. energy_pa: average atomic energy of the element's ground state structure after relaxation"
    )
    print(
        "3. a, c: lattice constants a and c of the element's ground state structure after relaxation"
    )
    print(
        "4. calculation details: psp=VASP.5.4.4 PBE recommended, ENCUT=400, KSPACING=0.15, EDIFF=1E-6, EDIFFG=-1E-2"
    )


parser = argparse.ArgumentParser(
    description="Common elements data.",
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

note()
