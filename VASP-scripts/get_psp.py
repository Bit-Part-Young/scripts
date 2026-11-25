#!/usr/bin/env python3

"""生成 VASP、pymatgen 推荐的 PBE 赝势 POTCAR 文件（不依赖 pymatgen）

赝势库根目录依次取: 环境变量 VASP_PSP_PATH、PMG_VASP_PSP_DIR、~/.pmgrc.yaml 中的
PMG_VASP_PSP_DIR。本文件中的赝势表和查找函数也被 psp_comparison.py 复用。
"""

import argparse
import gzip
import os
import re
import sys
from pathlib import Path

# Data source (pymatgen 已拆分，该文件现位于 pymatgen-core；2026-11-25 核对 v2026.10.31，仅 W 不同，见 PSP_NOTES):
# https://github.com/materialsproject/pymatgen-core/blob/main/src/pymatgen/io/vasp/MPRelaxSet.yaml
PBE_PYMATGEN = {
    "H": "H", "He": "He", "Li": "Li_sv", "Be": "Be_sv", "B": "B", "C": "C",
    "N": "N", "O": "O", "Mg": "Mg_pv", "Al": "Al", "Si": "Si", "P": "P",
    "S": "S", "Cl": "Cl", "Ca": "Ca_sv", "Sc": "Sc_sv", "Ti": "Ti_pv",
    "V": "V_pv", "Cr": "Cr_pv", "Mn": "Mn_pv", "Fe": "Fe_pv", "Co": "Co",
    "Ni": "Ni_pv", "Cu": "Cu_pv", "Zn": "Zn", "Zr": "Zr_sv", "Nb": "Nb_pv",
    "Mo": "Mo_pv", "Sn": "Sn_d", "Hf": "Hf_pv", "Ta": "Ta_pv", "W": "W_sv",
    "Pt": "Pt", "Au": "Au", "Pb": "Pb_d", "Ag": "Ag", "Pd": "Pd", "Na": "Na_pv",
    "U": "U",
}  # fmt: skip

# Data source: https://www.vasp.at/wiki/Recommended_pseudopotentials 中 potpaw.54, PBE（Standard）
# （2026-11-25 核对一致；potpaw.64 "before October 2026" 的推荐与此完全相同）
PBE_VASP = {
    "H": "H", "He": "He", "Li": "Li_sv", "Be": "Be", "B": "B", "C": "C",
    "N": "N", "O": "O", "Mg": "Mg", "Al": "Al", "Si": "Si", "P": "P",
    "S": "S", "Cl": "Cl", "Ca": "Ca_sv", "Sc": "Sc_sv", "Ti": "Ti_sv",
    "V": "V_sv", "Cr": "Cr_pv", "Mn": "Mn_pv", "Fe": "Fe", "Co": "Co",
    "Ni": "Ni", "Cu": "Cu", "Zn": "Zn", "Zr": "Zr_sv", "Nb": "Nb_sv",
    "Mo": "Mo_sv", "Sn": "Sn_d", "Hf": "Hf_pv", "Ta": "Ta_pv", "W": "W_sv",
    "Pt": "Pt", "Au": "Au", "Pb": "Pb_d", "Ag": "Ag", "Pd": "Pd", "Na": "Na_pv",
    "U": "U",
}  # fmt: skip

# 同上页面 potpaw.64, PBE 当前推荐（2026 年 10 月修订）：相对 .54 仅以下元素不同（2026-11-25 核对）
PBE_VASP64 = {
    **PBE_VASP,
    "Be": "Be_sv", "Mg": "Mg_pv", "Fe": "Fe_pv", "Co": "Co_pv", "Ni": "Ni_pv", "Cu": "Cu_pv",
}  # fmt: skip

PBE_RECOMMENDED = {"vasp": PBE_VASP, "vasp64": PBE_VASP64, "pymatgen": PBE_PYMATGEN}

# 与官方推荐不一致之处，运行时提示用户
PSP_NOTES = {
    ("pymatgen", "W"): "pymatgen MPRelaxSet recommends W_pv, which exists only in PBE.52 and older "
    "(not in PBE.54/64); W_sv is used instead.",
}

# 依次尝试的 PBE 子目录布局（本脚本原用 PAW_PBE，pymatgen 用 POT_GGA_PAW_PBE）
PBE_LAYOUTS = ("PAW_PBE/{}/POTCAR", "POT_GGA_PAW_PBE/POTCAR.{}", "potpaw_PBE/{}/POTCAR")


def psp_symbols(elements, recommended="vasp"):
    """元素符号 -> 推荐赝势符号；不支持的元素统一报错。"""
    table = PBE_RECOMMENDED[recommended]
    missing = [element for element in elements if element not in table]
    if missing:
        raise SystemExit(
            f"Error: no {recommended}-recommended PBE POTCAR for {', '.join(missing)}; "
            "please add it to PBE_VASP/PBE_PYMATGEN in get_psp.py."
        )
    for element in dict.fromkeys(elements):
        note = PSP_NOTES.get((recommended, element))
        if note:
            print(f"Note: {note}", file=sys.stderr)
    return [table[element] for element in elements]


def read_elements(structure_fn="POSCAR"):
    """读取 POSCAR 第 6 行元素符号（兼容 'Ti_sv/abc123' 这类带后缀的写法）。"""
    with open(structure_fn) as handle:
        lines = handle.readlines()
    tokens = lines[5].split() if len(lines) > 5 else []
    if not tokens or any(token[0].isdigit() for token in tokens):
        raise SystemExit(f"Error: {structure_fn} line 6 has no element symbols (VASP 4 format?).")
    return [token.split("/")[0].split("_")[0] for token in tokens]


def psp_root():
    """赝势库根目录。"""
    root = os.getenv("VASP_PSP_PATH") or os.getenv("PMG_VASP_PSP_DIR")
    pmgrc = Path.home() / ".pmgrc.yaml"
    if not root and pmgrc.is_file():
        match = re.search(r"^PMG_VASP_PSP_DIR:\s*(\S+)", pmgrc.read_text(), re.M)
        root = match and match.group(1).strip("'\"")
    if not root:
        raise SystemExit(
            "Error: set VASP_PSP_PATH (or PMG_VASP_PSP_DIR in ~/.pmgrc.yaml) to the POTCAR library root."
        )
    return Path(root).expanduser()


def find_potcar(psp_symbol, root=None):
    """在赝势库中查找单个元素 POTCAR（支持 .gz）。"""
    root = root or psp_root()
    for layout in PBE_LAYOUTS:
        for path in (root / layout.format(psp_symbol), root / (layout.format(psp_symbol) + ".gz")):
            if path.is_file():
                return path
    tried = ", ".join(layout.format(psp_symbol) for layout in PBE_LAYOUTS)
    raise SystemExit(f"Error: POTCAR for {psp_symbol} not found under {root} (tried {tried}).")


def read_potcar(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return handle.read()


def get_psp(structure_fn="POSCAR", psp_recommended="vasp", output="POTCAR"):
    """生成 VASP、pymatgen 推荐的 PBE 赝势 POTCAR 文件"""
    symbols = psp_symbols(read_elements(structure_fn), psp_recommended)
    root = psp_root()
    if psp_recommended == "vasp64":
        print(
            f"Note: vasp64 only picks the potpaw.64-recommended names; POTCAR files are read from {root}, "
            "make sure it is a potpaw.64 library if you need the .64 datasets.",
            file=sys.stderr,
        )
    # 先全部读取成功再写文件，避免中途失败留下残缺的 POTCAR
    text = "".join(read_potcar(find_potcar(symbol, root)) for symbol in symbols)
    Path(output).write_text(text)
    print(f"{output} ({' '.join(symbols)}) recommended by {psp_recommended} generated.")


def build_parser(description):
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument(
        "-pr",
        "--psp_recommended",
        default="vasp",
        choices=sorted(PBE_RECOMMENDED),
        help="recommended pseudopotential set: vasp = VASP potpaw.54 (= potpaw.64 before Oct 2026), "
        "vasp64 = VASP potpaw.64 (Oct 2026 revision), "
        "pymatgen = MPRelaxSet",
    )
    parser.add_argument("structure_fn", nargs="?", default="POSCAR", help="structure filename")
    return parser


if __name__ == "__main__":
    args = build_parser(
        "Generate VASP, pymatgen recommended PBE pseudopotentials (pymatgen independent)."
    ).parse_args()
    get_psp(args.structure_fn, args.psp_recommended)
