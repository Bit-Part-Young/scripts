"""固定（特定）原子层 x/y/z 轴"""

import numpy as np
import pandas as pd
from pymatgen.core.structure import Structure


def fix_layer(
    structure: Structure,
    layer_index: list[int] | None = None,
    fixed_axes: list[str] = ["T", "T", "T"],
    precision: float = 0.001,
) -> Structure:
    """
    固定（特定）原子层 x/y/z 轴
    """

    natoms = structure.num_sites
    positions = structure.cart_coords
    # 移除原子坐标轴原有的固定信息
    # structure.remove_site_property(property_name="selective_dynamics")
    # site_properties = structure.site_properties

    z_coords = positions[:, 2]
    z_coords_rounded = precision * np.round(z_coords / precision)

    z_unique = np.unique(z_coords_rounded)

    print(f"Number of atom layer: {len(z_unique)}.\n")

    # 获取每个原子所在的原子层
    layer_index_list = []
    for z_coord in z_coords_rounded:
        for axis_index, z in enumerate(z_unique, start=1):
            if np.isclose(z, z_coord):
                layer_index_list.append(axis_index)

    data = {
        "x": positions[:, 0],
        "y": positions[:, 1],
        "z": positions[:, 2],
        "layer_index": layer_index_list,
    }
    df = pd.DataFrame(data)

    if layer_index is None:
        fixed_layer_index = df.index

        print(f"All {len(fixed_layer_index)} atoms are fixed.")
    elif isinstance(layer_index, list):
        fixed_layer_index = df[df["layer_index"].isin(layer_index)].index

        print(
            f"No. {layer_index} layer(s), total {len(fixed_layer_index)} atoms are fixed."
        )

    selective_dynamics = np.ones((natoms, 3))

    for axis_index, fixed_string in enumerate(fixed_axes):
        if fixed_string == "F":
            selective_dynamics[fixed_layer_index, axis_index] = 0

    print(f"Axes fixed: {fixed_axes}.\n")

    structure.add_site_property("selective_dynamics", selective_dynamics)

    return structure


if __name__ == "__main__":

    # 示例
    structure_fn = "POSCAR"
    structure = Structure.from_file(structure_fn)

    # 移除原子坐标轴原有的固定信息
    structure.remove_site_property(property_name="selective_dynamics")

    # 固定所有原子层
    layer_index = None
    # 固定单个原子层
    # layer_index = [1]
    # 固定多个原子层
    # layer_index = [1, 2]
    # 固定 x 轴
    fixed_axis = ["F", "T", "T"]
    # 固定 x、y 轴
    # fixed_axis = ["F", "F", "T"]

    structure_fixed = fix_layer(
        structure=structure,
        layer_index=layer_index,
        fixed_axes=fixed_axis,
    )

    structure_fixed.to("POSCAR_fixed", fmt="poscar")
