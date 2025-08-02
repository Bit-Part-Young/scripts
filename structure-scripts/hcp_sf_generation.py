#!/usr/bin/env python3

"""生成 HCP I1, I2, T2, E 型 4 种层错构型"""

import argparse
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.io.vasp import Poscar
from ase.io import write
from ase.atoms import Atoms
import numpy as np


def save_poscar(
    structure: Structure | Atoms,
    output_fn: str,
    sort: bool = True,
    direct: bool = True,
):
    """Save structure to VASP POSCAR file

    Args:
        structure: pymatgen Structure or ASE Atoms object
        sort: whether to sort the structure
        direct: whether to use direct coordinates
        output_fn: output filename
    """

    if isinstance(structure, Structure):
        poscar = Poscar(structure, sort_structure=sort)
        poscar.write_file(output_fn, direct=direct, significant_figures=10)
    elif isinstance(structure, Atoms):
        write(output_fn, images=structure, format="vasp", direct=direct, sort=sort)


def wrap_pos(structure: Structure) -> Structure:
    """Wrap the position of the structure to 0-1.

    Args:
        structure: pymatgen Structure object
    """

    atoms = structure.to_ase_atoms()
    atoms.wrap()

    structure = Structure.from_ase_atoms(atoms)

    return structure


def add_vacuum(
    structure: Structure | str,
    vacuum: float,
    mode: str,
    to_poscar: bool = False,
    to_print: bool = False,
) -> Structure:
    """Add vacuums to top, bottom, or both side in z direction

    Args:
        structure: Structure object or structure file
        vacuum: vacuum thickness (Å)
        mode: vacuum mode, top, bottom, both
        to_poscar: whether to save to file
        to_print: whether to print the structure

    Returns:
        Structure object
    """

    if isinstance(structure, str):
        structure = Structure.from_file(structure)

    species = structure.species
    lattice = structure.lattice.matrix
    lattice_copy = lattice.copy()
    positions = structure.cart_coords

    if mode == "top":
        lattice_copy[2, 2] += vacuum
    elif mode == "bottom":
        lattice_copy[2, 2] += vacuum
        positions[:, 2] += vacuum
    elif mode == "both":
        lattice_copy[2, 2] += vacuum * 2
        positions[:, 2] += vacuum

    structure_vacuum = Structure(
        lattice=lattice_copy,
        species=species,
        coords=positions,
        coords_are_cartesian=True,
    )

    if to_print:
        print(f"\nAdd {vacuum} Å vacuum to {mode} side.\n")
        print(structure_vacuum)

    if to_poscar:
        output_fn = f"{structure.composition.reduced_formula}_vacuum_{mode}.vasp"
        save_poscar(structure_vacuum, output_fn)

        print(f"\n{output_fn} saved.")

    return structure_vacuum


def hcp_orthogonal_unitcell(element: str, a: float, c: float) -> Structure:
    """Generate HCP orthogonal unit cell

    Args:
        element: element symbol
        a: lattice constant a
        c: lattice constant c

    Returns:
        Structure object
    """

    # x 方向为 <10-10>
    structure = Structure(
        lattice=Lattice.from_parameters(
            a=1.7320508076 * a, b=a, c=c, alpha=90.0, beta=90.0, gamma=90.0
        ),
        species=[element] * 4,
        coords=[
            [0.50000000, 0.50000000, 0.00000000],
            [0.00000000, 0.00000000, 0.00000000],
            [0.66666667, 0.00000000, 0.50000000],
            [0.16666667, 0.50000000, 0.50000000],
        ],
    )

    return structure


def sf_I2(
    structure: Structure,
    duplicate: list[int] = [1, 2, 6],
    vacuum: float = 0.0,
    to_poscar: bool = True,
) -> Structure:
    """生成 I2 型层错构型"""

    structure_supercell: Structure = structure * duplicate
    structure_copy = structure_supercell.copy()

    burger_vector_mod = structure.lattice.a

    # 上半部分移动 1/3 <10-10>
    z_half = structure_supercell.lattice.c * 0.48
    vector = np.array([burger_vector_mod / 3, 0, 0])
    indices = np.where(structure_copy.cart_coords[:, 2] > z_half)[0]
    structure_copy.translate_sites(
        indices, vector, frac_coords=False, to_unit_cell=False
    )

    structure_copy = wrap_pos(structure_copy)
    structure_copy = add_vacuum(structure_copy, vacuum=vacuum, mode="top")

    site_property = np.zeros((structure_copy.num_sites, 3), dtype=int)
    site_property[:, -1] = 1
    structure_copy.add_site_property("selective_dynamics", site_property.tolist())

    if to_poscar:
        structure_fn = "sf_I2.vasp"
        save_poscar(structure_copy, structure_fn)
        print(f"\nHCP I2 type stacking fault configuration {structure_fn} saved.")

    return structure_copy


def sf_T2(
    structure: Structure,
    duplicate: list[int] = [1, 2, 6],
    vacuum: float = 0.0,
    to_poscar: bool = True,
) -> Structure:
    """生成 T2 型层错构型"""

    structure_supercell: Structure = sf_I2(structure, duplicate, to_poscar=False)
    structure_copy = structure_supercell.copy()

    burger_vector_mod = structure.lattice.a

    # 在 I2 型层错的基础上，C 以上的所有原子层移动 -1/3 <10-10>
    z_half = structure_supercell.lattice.c * 0.52
    vector = np.array([-burger_vector_mod / 3, 0, 0])
    indices = np.where(structure_copy.cart_coords[:, 2] > z_half)[0]
    structure_copy.translate_sites(
        indices, vector, frac_coords=False, to_unit_cell=False
    )

    structure_copy = wrap_pos(structure_copy)
    structure_copy = add_vacuum(structure_copy, vacuum=vacuum, mode="top")

    site_property = np.zeros((structure_copy.num_sites, 3), dtype=int)
    site_property[:, -1] = 1
    structure_copy.add_site_property("selective_dynamics", site_property.tolist())

    if to_poscar:
        structure_fn = "sf_T2.vasp"
        save_poscar(structure_copy, structure_fn)
        print(f"\nHCP T2 type stacking fault configuration {structure_fn} saved.")

    return structure_copy


def sf_E(
    structure: Structure,
    duplicate: list[int] = [1, 2, 6],
    vacuum: float = 0.0,
    to_poscar: bool = True,
) -> Structure:
    """生成 E 型层错构型"""

    # 元素、层间距、伯氏矢量
    specie = structure.sites[0].species_string
    layer_distance = structure.lattice.c / 2
    burger_vector_mod = structure.lattice.a

    structure_supercell: Structure = structure * duplicate
    structure_copy = structure_supercell.copy()
    natoms = structure_supercell.num_sites
    a, b, c = structure_supercell.lattice.abc
    positions_cartesian = structure_copy.cart_coords

    # 将上半部分的原子 z 坐标增加一个层间距
    z_half = structure_supercell.lattice.c * 0.48
    indices = np.where(structure_copy.cart_coords[:, 2] > z_half)[0]
    positions_cartesian[indices, 2] += layer_distance

    # 获取第 2 个原子层的 B 类堆垛原子 x、y 坐标，使其移动 1/6 <10-10>，形成 C 类堆垛原子并置于胞的中心（CNA 识别中间层只有 1 层为 FCC 类原子）
    # 获取第 1 个原子层的 A 类堆垛原子 x、y 坐标，使其移动 1/3 <10-10>，形成 C 类堆垛原子并置于胞的中心（CNA 识别中间层有 3 层为 FCC 类原子）
    # [ ] 哪个对？目前采用的第一种
    indices2 = np.where(
        np.isclose(positions_cartesian[:, 2], layer_distance, atol=1e-6)
    )[0]
    position_cartesian_new = positions_cartesian[indices2]
    position_cartesian_new += np.array(
        [burger_vector_mod / 6, 0, c / 2 - layer_distance]
    )
    natoms_new = position_cartesian_new.shape[0]

    positions_cartesian = np.concatenate([positions_cartesian, position_cartesian_new])

    # c 增加一个层间距
    structure_new = Structure(
        lattice=Lattice.from_parameters(
            a=a, b=b, c=c + layer_distance, alpha=90.0, beta=90.0, gamma=90.0
        ),
        species=[specie] * (natoms + natoms_new),
        coords=positions_cartesian,
        coords_are_cartesian=True,
    )

    structure_new = wrap_pos(structure_new)
    structure_new = add_vacuum(structure_new, vacuum=vacuum, mode="top")

    site_property = np.zeros((structure_new.num_sites, 3), dtype=int)
    site_property[:, -1] = 1
    structure_new.add_site_property("selective_dynamics", site_property.tolist())

    if to_poscar:
        structure_fn = "sf_E.vasp"
        save_poscar(structure_new, structure_fn)
        print(f"\nHCP E type stacking fault configuration {structure_fn} saved.")

    return structure_new


def sf_I1(
    structure: Structure,
    duplicate: list[int] = [1, 2, 6],
    vacuum: float = 0.0,
    to_poscar: bool = True,
) -> Structure:
    """生成 I1 型层错构型"""

    # 层间距、伯氏矢量
    layer_distance = structure.lattice.c / 2
    burger_vector_mod = structure.lattice.a

    structure_supercell: Structure = structure * duplicate
    structure_copy = structure_supercell.copy()
    a, b, c = structure_supercell.lattice.abc
    positions_cartesian = structure_copy.cart_coords

    # 删除上半部分的第一层原子
    indices2 = np.where(np.isclose(positions_cartesian[:, 2], c / 2, atol=1e-6))[0]
    structure_copy.remove_sites(indices2)

    # 移动上半部分的原子并使其 z 坐标减少一个层间距
    z_half = structure_supercell.lattice.c * 0.48
    vector = np.array([burger_vector_mod / 6, 0, -layer_distance])
    indices = np.where(structure_copy.cart_coords[:, 2] > z_half)[0]
    structure_copy.translate_sites(
        indices, vector, frac_coords=False, to_unit_cell=False
    )

    # c 减去一个层间距
    structure_new = Structure(
        lattice=Lattice.from_parameters(
            a=a, b=b, c=c - layer_distance, alpha=90.0, beta=90.0, gamma=90.0
        ),
        species=structure_copy.species,
        coords=structure_copy.cart_coords,
        coords_are_cartesian=True,
    )

    structure_new = wrap_pos(structure_new)
    structure_new = add_vacuum(structure_new, vacuum=vacuum, mode="top")

    site_property = np.zeros((structure_new.num_sites, 3), dtype=int)
    site_property[:, -1] = 1
    structure_new.add_site_property("selective_dynamics", site_property.tolist())

    if to_poscar:
        structure_fn = "sf_I1.vasp"
        save_poscar(structure_new, structure_fn)
        print(f"\nHCP I1 type stacking fault configuration {structure_fn} saved.")

    return structure_new


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate HCP I1, I2, T2, E four types stacking fault configurations.",
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "-e",
        "--element",
        metavar="element",
        help="element symbol",
    )
    parser.add_argument(
        "-lc",
        "--lattice_constants",
        type=float,
        nargs=2,
        metavar="lattice_constants",
        help="lattice constants",
    )
    parser.add_argument(
        "-t",
        "--type",
        metavar="type",
        choices=["I1", "I2", "T2", "E"],
        help="stacking fault type (I1, I2, T2, E)",
    )
    parser.add_argument(
        "-d",
        "--duplicate",
        nargs=3,
        type=int,
        default=[1, 2, 6],
        metavar="duplicate",
        help="duplicate structure in the three directions (default: 1 2 6)",
    )
    parser.add_argument(
        "-vac",
        "--vacuum",
        type=float,
        default=0.0,
        metavar="vacuum",
        help="vacuum thickness (default: 0.0 Å)",
    )
    parser.add_argument(
        "-o",
        action="store_true",
        default=False,
        help="save to POSCAR file",
    )

    args = parser.parse_args()

    structure = hcp_orthogonal_unitcell(
        args.element, args.lattice_constants[0], args.lattice_constants[1]
    )

    if args.type == "I2":
        sf_I2(structure, args.duplicate, args.vacuum, args.o)
    elif args.type == "T2":
        sf_T2(structure, args.duplicate, args.vacuum, args.o)
    elif args.type == "E":
        sf_E(structure, args.duplicate, args.vacuum, args.o)
    elif args.type == "I1":
        sf_I1(structure, args.duplicate, args.vacuum, args.o)
    else:
        raise ValueError(f"Invalid stacking fault type: {args.type}")
