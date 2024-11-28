"""固定（特定）原子层 x/y/z 轴"""

import numpy as np
import pandas as pd
from pymatgen.core.structure import Structure


def fix_layer(
    structure: Structure,
    layer_index: int | list[int] | None = None,
    fixed_axis: int | list[int] | None = None,
    precision: float = 0.001,
    to_poscar: bool = False,
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
        for i, z in enumerate(z_unique, start=1):
            if np.isclose(z, z_coord):
                layer_index_list.append(i)

    data = {
        "x": positions[:, 0],
        "y": positions[:, 1],
        "z": positions[:, 2],
        "layer_index": layer_index_list,
    }
    df = pd.DataFrame(data)

    if layer_index is None:
        fixed_layer_index = df.index

        print("All atoms are fixed.")
    elif isinstance(layer_index, int):
        fixed_layer_index = df[df["layer_index"] == layer_index].index

        print(
            f"No. {layer_index} layer, total {len(fixed_layer_index)} atoms are fixed."
        )
    elif isinstance(layer_index, list):
        fixed_layer_index = df[df["layer_index"].isin(layer_index)].index

        print(
            f"No. {layer_index} layers, total {len(fixed_layer_index)} atoms are fixed."
        )

    selective_dynamics = np.ones((natoms, 3))

    if isinstance(fixed_axis, int):
        selective_dynamics[fixed_layer_index, fixed_axis] = 0
    elif isinstance(fixed_axis, list):
        selective_dynamics[np.ix_(fixed_layer_index, fixed_axis)] = 0

    structure.add_site_property("selective_dynamics", selective_dynamics)

    if to_poscar:
        structure.to("POSCAR", fmt="poscar")

    return structure


if __name__ == "__main__":

    # 示例
    structure_fn = "POSCAR_1"
    structure = Structure.from_file(structure_fn)

    # 移除原子坐标轴原有的固定信息
    structure.remove_site_property(property_name="selective_dynamics")

    # 固定单个原子层
    layer_index = 1
    # 固定多个原子层
    layer_index = [1, 2]
    # 固定 x 轴
    fixed_axis = 0
    # 固定 x、y 轴
    fixed_axis = [0, 1]

    fix_layer(
        structure=structure,
        layer_index=layer_index,
        fixed_axis=fixed_axis,
        to_poscar=True,
    )
