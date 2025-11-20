#!/usr/bin/env python3

"""获取 Materials Project 中的一元、二元数据"""

import argparse
import os

import pandas as pd
from mp_api.client import MPRester
from pymatgen.core.structure import Structure

# 注: energy_per_atom 的数值非 VASP 计算值
fields = [
    "material_id",
    "structure",
    "composition_reduced",
    "formula_pretty",
    "symmetry",
    "nsites",
    "energy_per_atom",
    "formation_energy_per_atom",
    "energy_above_hull",
    "is_stable",
    "structure",
]

API_KEY = os.getenv("PMG_MAPI_KEY")
if API_KEY is None:
    raise OSError("PMG_MAPI_KEY environment variable is not set. Please check!")

# 优先以元素 Al 作为 0K 二元相图的 x 轴变量
element_sequence_list = ["Al", "Ti", "Nb", "Mo", "Zr", "V"]
# element_sequence_list = ["Ti", "Al", "Nb", "Mo", "Zr", "V"]


def get_mp_data(
    elements: str | list[str],
    stable: bool | None = None,
    energy_above_hull: float | None = None,
):
    """获取 Materials Project 中的一元、二元数据"""

    if isinstance(elements, list):
        # eg. ["Ti", "Al", "Ti-Al"]
        chemsys_list = elements + ["-".join(elements)]
    else:
        chemsys_list = elements

    if energy_above_hull is not None:
        energy_above_hull = (0.0, energy_above_hull)

    if stable is False:
        stable = None

    with MPRester(API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            chemsys=chemsys_list,
            fields=fields,
            is_stable=stable,
            energy_above_hull=energy_above_hull,
        )

    ndocs = len(docs)

    if isinstance(elements, list):
        print(f"\nTotally get {ndocs} {' '.join(elements)} MP binary data.\n")
    else:
        print(f"\nTotally get {ndocs} {elements} MP unary data.\n")

    data_list = []
    for doc in docs:

        structure: Structure = doc.structure
        # 删除磁性信息
        if "magmom" in structure.site_properties.keys():
            structure.remove_site_property("magmom")

        atoms = structure.to_ase_atoms()
        natoms = len(atoms)
        formula = atoms.get_chemical_formula()
        cell_info = atoms.cell.cellpar().round(3).tolist()

        # formula 直接使用 doc.formula_pretty 会导致后面的 0K 二元相图绘制不准确
        # formula_pretty 会导致获取的所有一元的化学式均为 如 Ti Al 之类，而非 Ti2 Al4
        # formula_pretty 对金属间化合物无影响
        data_dict = {
            "formula": formula,
            "natoms": natoms,
            "material_id": str(doc.material_id),
            "energy": doc.energy_per_atom * doc.nsites,
            "fe": doc.formation_energy_per_atom,
            "e_above_hull": doc.energy_above_hull,
            "is_stable": doc.is_stable,
            "crystal_system": str(doc.symmetry.crystal_system),
            "spg_symbol": doc.symmetry.symbol,
            "spg_number": doc.symmetry.number,
            "cell_info": cell_info,
        }

        if isinstance(elements, list):
            # 按照给定的元素顺序排序
            elements.sort(key=lambda x: element_sequence_list.index(x))

            composition: dict[str, float] = doc.composition_reduced.as_dict()

            # 成分 分数形式
            composition_fractional = {
                element: composition.get(element, 0.0) / sum(composition.values())
                for element in elements
            }

            data_dict.update(composition_fractional)

        data_list.append(data_dict)

    df = pd.DataFrame(data_list).round(5)
    print()
    print(df)

    if isinstance(elements, list):
        csv_fn = f"{'_'.join(elements)}_mp.csv"
    else:
        csv_fn = f"{elements}_mp.csv"
    df.to_csv(csv_fn, index=False)

    print("\nNote: The energy in MP is not real DFT calculated energy?")
    print("Note: The stable entry in MP for pure element may be not the common one.")

    print(f"\nData saved to {csv_fn}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get unary/binary data from Materials Project.",
        epilog="Author: SLY.",
    )

    parser.add_argument("elements", nargs="+", help="elements (e.g. Ti, Ti Al)")
    parser.add_argument("--stable", action="store_true", help="only get stable data")
    parser.add_argument(
        "-eah",
        "--energy_above_hull",
        type=float,
        metavar="FLOAT",
        help="maximum energy above hull",
    )

    args = parser.parse_args()

    if len(args.elements) == 1:
        args.elements = args.elements[0]
    elif len(args.elements) > 2:
        raise ValueError("Number of elements should not more than 2.")

    get_mp_data(
        elements=args.elements,
        stable=args.stable,
        energy_above_hull=args.energy_above_hull,
    )
