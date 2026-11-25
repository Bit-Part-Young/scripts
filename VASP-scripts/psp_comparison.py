#!/usr/bin/env python3

"""比较 VASP(potpaw.54/.64) 和 pymatgen 推荐的 PBE 赝势（TITEL、ENMAX、ZVAL、VRHFIN），仅依赖标准库"""

import argparse
import re

from get_psp import find_potcar, psp_root, psp_symbols, read_potcar

COLUMNS = ("", "POTCAR", "TITEL", "ENMAX", "ZVAL", "VRHFIN")


def get_potcar_info(psp_symbol, root=None):
    """读取单元素 POTCAR 头部信息"""
    text = read_potcar(find_potcar(psp_symbol, root))

    def field(pattern):
        return re.search(pattern, text, re.M).group(1).strip()

    return {
        "POTCAR": psp_symbol,
        "TITEL": field(r"^\s*TITEL\s*=\s*(.+)$"),
        "ENMAX": "%.3f" % float(field(r"ENMAX\s*=\s*([-\d.]+)")),
        "ZVAL": "%g" % float(field(r"ZVAL\s*=\s*([-\d.]+)")),
        # 部分官方 POTCAR（如 Ni、Pb、Bi）的 VRHFIN 本身就为空
        "VRHFIN": field(r"VRHFIN\s*=\s*[^:]*:(.*)$") or "-",
    }


def potcar_pbe_compare(element_symbols):
    """比较 VASP(.54/.64) 和 pymatgen 推荐的 PBE 赝势；推荐相同的合并为一行"""
    root = psp_root()
    sources = {"VASP.54": "vasp", "VASP.64": "vasp64", "pymatgen": "pymatgen"}
    symbols = {label: psp_symbols(element_symbols, key) for label, key in sources.items()}
    rows = []
    for i, element in enumerate(element_symbols):
        groups = {}
        for label in sources:
            groups.setdefault(symbols[label][i], []).append(label)
        for psp, labels in groups.items():
            rows.append({"": f"{element} {'='.join(labels)}", **get_potcar_info(psp, root)})

    widths = {c: max(len(c), *(len(row[c]) for row in rows)) for c in COLUMNS}
    numeric = {"ENMAX", "ZVAL"}
    for row in [dict(zip(COLUMNS, COLUMNS)), *rows]:
        print("  ".join(row[c].rjust(widths[c]) if c in numeric else row[c].ljust(widths[c])
                        for c in COLUMNS).rstrip())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare VASP (potpaw.54/.64) and pymatgen recommended PBE pseudopotentials.",
        epilog="Example: %(prog)s Ti Al W. Author: SLY.",
    )
    parser.add_argument("element_symbols", nargs="+", help="element symbols, eg. Ti Al")
    potcar_pbe_compare(parser.parse_args().element_symbols)
