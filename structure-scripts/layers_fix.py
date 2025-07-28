"""固定（特定）原子层 x/y/z 轴"""

import argparse

import numpy as np
import pandas as pd
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar


def layers_fix(
    structure_fn: str,
    layer_indices: int | list[int] | None = None,
    fixed_axes: list[str] = ["T", "T", "T"],
    to_poscar: bool = False,
    precision: float = 0.001,
) -> Structure:
    """固定（特定）原子层 x/y/z 轴"""

    structure = Structure.from_file(structure_fn)

    natoms = structure.num_sites
    positions = structure.cart_coords

    positions_z = positions[:, 2]
    positions_z_rounded = precision * np.round(positions_z / precision)
    positions_z_unique = np.unique(positions_z_rounded)

    print(f"\nAtoms count: {natoms}; Atomic layers count: {len(positions_z_unique)}.\n")

    # 获取每个原子所在的原子层索引
    layer_index_list = []
    for z_round in positions_z_rounded:
        for i, z_unique in enumerate(positions_z_unique, start=1):
            if np.isclose(z_unique, z_round):
                layer_index_list.append(i)

    data = {
        "x": positions[:, 0],
        "y": positions[:, 1],
        "z": positions[:, 2],
        "layer_index": layer_index_list,
    }
    df = pd.DataFrame(data)

    if layer_indices is None:
        fixed_layer_indices = df.index

        print(f"All {len(fixed_layer_indices)} atoms are fixed.")
    elif isinstance(layer_indices, list):
        fixed_layer_indices = df[df["layer_index"].isin(layer_indices)].index

        print(
            f"No. {layer_indices} layer(s), total {len(fixed_layer_indices)} atoms are fixed."
        )

    selective_dynamics = structure.site_properties.get("selective_dynamics", None)
    if selective_dynamics is None:
        selective_dynamics = np.ones((natoms, 3), dtype=int)
    else:
        selective_dynamics = np.array(selective_dynamics)

    for axis_index, fixed_string in enumerate(fixed_axes):
        if fixed_string == "F":
            selective_dynamics[fixed_layer_indices, axis_index] = 0

    print(f"Axes fixed: {fixed_axes}.")

    structure.add_site_property("selective_dynamics", selective_dynamics.tolist())

    if to_poscar:
        output_fn = "layers_fixed.vasp"
        poscar = Poscar(structure, sort_structure=True)
        poscar.write_file(output_fn, direct=True, significant_figures=10)

        print(f"\n{output_fn} saved.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fix atomic layers in x, y, z axis",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
    )
    parser.add_argument(
        "structure_fn",
        nargs="?",
        default="POSCAR",
        metavar="structure_fn",
        help="structure filename",
    )
    parser.add_argument(
        "-li",
        "--layer_indices",
        type=int,
        nargs="*",
        metavar="layer_indices",
        help="atomic layer indices (eg. 1, 1 2)",
    )
    parser.add_argument(
        "-fa",
        "--fixed_axes",
        nargs=3,
        metavar="fixed_axes",
        help="fixed axes (eg. F F T, F F F)",
    )
    parser.add_argument("-o", action="store_true", help="write file")

    args = parser.parse_args()

    layers_fix(
        structure_fn=args.structure_fn,
        layer_indices=args.layer_indices,
        fixed_axes=args.fixed_axes,
        to_poscar=args.o,
    )
