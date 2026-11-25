#!/usr/bin/env python3

"""
元素周期表常见元素数据

Z: 原子序数
mass: 相对原子质量
r(Å): 原子半径（无统一定值，不一定准确，取自 pymatgen）
T_m(K): 熔点
ρ(g/cm^3): 密度
E_pa(eV/atom): 基态结构弛豫后的平均原子能量
structure: 基态结构
a(Å), c(Å): 基态结构弛豫后的晶格常数
mp-id: Materials Project ID
"""

import argparse
import sys

import pandas as pd

COLUMNS = ["Z", "mass", "r(Å)", "T_m(K)", "ρ(g/cm^3)", "E_pa(eV/atom)", "structure", "a(Å)", "c(Å)", "mp-id"]

DATA = {
    "Ti": (22, 47.867, 1.4, 1941.0, 4.51, -7.8, "hcp", 2.93, 4.65, "mp-46"),
    "Al": (13, 26.982, 1.25, 933.47, 2.70, -3.7, "fcc", 4.04, 4.04, "mp-134"),
    "Nb": (41, 92.906, 1.45, 2750.0, 8.57, -10.2, "bcc", 3.31, 3.31, "mp-75"),
    "Zr": (40, 91.224, 1.55, 2128.0, 6.51, -8.5, "hcp", 3.23, 5.15, "mp-131"),
    "V": (23, 50.941, 1.35, 2183.0, 6.11, -8.9, "bcc", 2.98, 2.98, "mp-146"),
    "Mo": (42, 95.94, 1.45, 2896.0, 10.28, -10.9, "bcc", 3.16, 3.16, "mp-129"),
    "Ta": (73, 180.948, 1.45, 3290.0, 16.65, -11.8, "bcc", 3.32, 3.32, "mp-50"),
    "W": (74, 183.84, 1.35, 3695.0, 19.25, -13.0, "bcc", 3.18, 3.18, "mp-91"),
    "Hf": (72, 178.49, 1.55, 2506.0, 13.31, -9.9, "hcp", 3.20, 5.06, "mp-103"),
    "Sn": (50, 118.71, 1.45, 505.08, 7.31, -3.8, "diamond", 6.65, 6.65, "mp-117"),
    "Mg": (12, 24.305, 1.5, 923.0, 1.74, -1.5, "hcp", 3.20, 5.16, "mp-153"),
    "Cu": (29, 63.546, 1.35, 1357.77, 8.92, -3.7, "fcc", 3.63, 3.63, "mp-30"),
    "Si": (14, 28.0855, 1.1, 1687.0, 2.33, -5.4, "diamond", 5.47, 5.47, "mp-149"),
    "Sc": (21, 44.956, 1.6, 1814.0, 2.99, -6.3, "hcp", 3.32, 5.16, "mp-67"),
    "Ni": (28, 58.693, 1.35, 1728.0, 8.91, -5.5, "fcc", 3.52, 3.52, "mp-23"),
    "Cr": (24, 51.996, 1.4, 2180.0, 7.14, -9.5, "bcc", 2.84, 2.84, "mp-90"),
    "Ag": (47, 107.868, 1.6, 1234.93, 10.49, -2.7, "fcc", 4.14, 4.14, "mp-124"),
    "Au": (79, 196.967, 1.35, 1337.33, 19.30, -3.2, "fcc", 4.16, 4.16, "mp-81"),
    "Pb": (82, 207.2, 1.8, 600.61, 11.34, None, None, None, None, "mp-20483"),
    "Pd": (46, 106.42, 1.4, 1828.05, 12.02, -5.2, "fcc", 3.94, 3.94, "mp-2"),
    "Pt": (78, 195.084, 1.35, 2041.4, 21.09, -6.1, "fcc", 3.97, 3.97, "mp-126"),
}

NOTE = """
Note:
1. r: atomic radius (no unified value, not necessarily accurate), from pymatgen.
2. E_pa: average energy per atom of the relaxed ground-state structure.
3. a, c: lattice constants of the relaxed ground-state structure.
4. Calculation details: VASP.5.4.4 PBE recommended PSPs, ENCUT=400, KSPACING=0.15,
   EDIFF=1E-6, EDIFFG=-1E-2.
"""


def element_data(elements: list[str] | None = None, sort: bool = False) -> pd.DataFrame:
    """返回元素数据表；elements 为 None 时返回全部元素"""

    df = pd.DataFrame.from_dict(DATA, orient="index", columns=COLUMNS)

    if elements:
        elements = [e.capitalize() for e in elements]
        missing = [e for e in elements if e not in df.index]
        if missing:
            sys.exit(f"Element(s) {missing} not found. Available: {' '.join(df.index)}")
        df = df.loc[elements]

    return df.sort_values("Z") if sort else df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Show data of common elements.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument(
        "-e", "--element", nargs="+", metavar="EL", help="element symbols (e.g. Ti Al); omit to show all"
    )
    parser.add_argument("--sort", action="store_true", help="sort by atomic number")
    args = parser.parse_args()

    print()
    print(element_data(args.element, args.sort).to_string(na_rep="-"))
    print(NOTE)
