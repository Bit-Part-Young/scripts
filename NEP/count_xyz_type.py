#!/usr/bin/env python3

"""统计 xyz 文件中的元素类别及其数量"""

import argparse

elements_reference = ["Ti", "Al", "Nb", "Mo", "Zr", "V"]


def count_xyz_type(filename: str) -> dict:
    """统计 xyz 文件中的元素类别及其数量"""

    with open(filename, "r") as file:
        lines = file.readlines()

    unique_elements_list = []
    i = 0
    while i < len(lines):
        natoms = int(lines[i].strip())
        atoms = lines[i + 2 : i + 2 + natoms]

        elements = [atom.split()[0] for atom in atoms]
        unique_elements = list(set(elements))
        # 按照 elements_reference 的顺序排序
        unique_elements.sort(key=lambda x: elements_reference.index(x))
        elements_str = "-".join(unique_elements)

        unique_elements_list.append(elements_str)

        i += 2 + natoms

    # 统计 elements_str 的类型及数量
    count_dict = {}
    for elements_str in unique_elements_list:
        if elements_str in count_dict:
            count_dict[elements_str] += 1
        else:
            count_dict[elements_str] = 1

    return count_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count element type in xyz file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY",
    )

    parser.add_argument("xyz_fn", type=str, help="xyz filename")

    args = parser.parse_args()
    xyz_fn = args.xyz_fn

    count_dict = count_xyz_type(xyz_fn)

    print(f"\n{xyz_fn} count info:")
    for elements_str, count in count_dict.items():
        print(f"{elements_str}: {count}")
    print(f"Total: {sum(count_dict.values())}")
