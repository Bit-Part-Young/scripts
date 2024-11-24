"""固定（特定）原子层 z 轴"""

import numpy as np
import pandas as pd
from pymatgen.core.structure import Structure


def fix_sf_layer(
    structure: Structure,
    layer_number: int | list[int] | None = None,
    precision: float = 0.001,
    to_poscar: bool = False,
) -> Structure:
    """固定（特定）原子层 z 轴"""

    natoms = structure.num_sites
    positions = structure.cart_coords
    # 移除原子坐标轴原有的固定信息
    # structure.remove_site_property(property_name="selective_dynamics")
    # site_properties = structure.site_properties

    z_coords = positions[:, 2]
    z_coords_rounded = precision * np.round(z_coords / precision)

    z_unique = np.unique(z_coords_rounded)
    z_unique_sorted = np.sort(z_unique)

    print(f"Number of atom layer: {len(z_unique)}.\n")

    # 获取每个原子所在的原子层
    layer_num_list = []
    for z_coord in z_coords_rounded:
        for i, z in enumerate(z_unique_sorted, start=1):
            if np.isclose(z, z_coord):
                layer_num_list.append(i)

    data = {
        "x": positions[:, 0],
        "y": positions[:, 1],
        "z": positions[:, 2],
        "layer_number": layer_num_list,
    }
    df = pd.DataFrame(data)

    if layer_number is None:
        fix_index = df.index
        print("z axis of all layers are fixed.")
    elif isinstance(layer_number, int):
        fix_index = df[df["layer_number"] == layer_number].index
        print(f"z axis of No. {fix_index} layer is fixed.")
    elif isinstance(layer_number, list):
        fix_index = df[df["layer_number"].isin(layer_number)].index
        print(f"z axis of No. {fix_index} layers are fixed.")

    selective_dynamics = np.ones((natoms, 3))
    # 固定选中原子层的 z 轴
    selective_dynamics[fix_index, 2] = 0

    structure.add_site_property("selective_dynamics", selective_dynamics)

    if to_poscar:
        structure.to("POSCAR", fmt="poscar")

    return structure


if __name__ == "__main__":

    structure_fn = "POSCAR_1"
    structure = Structure.from_file(structure_fn)
    layer_number = 6
    layer_number = [6, 7]

    fix_sf_layer(
        structure=structure,
        layer_number=layer_number,
        to_poscar=True,
    )
