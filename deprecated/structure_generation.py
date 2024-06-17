#!/usr/bin/env python3

import argparse
import math

from ase.io.lammpsdata import write_lammps_data
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def supercell_gen(struc, x_scale, y_scale, z_scale, struc_type):
    # sp = SpacegroupAnalyzer(struc)
    struc_copy = struc.copy()
    struc_copy.make_supercell([x_scale, y_scale, z_scale])
    struc_ase = AseAtomsAdaptor.get_atoms(struc_copy)

    pos = Poscar(struc_copy)
    pos.write_file("POSCAR")

    write_lammps_data("data.pos", struc_ase, atom_style="atomic")

    print(f"{struc_type} structure vasp POSCAR and lammps pos file are generated.\n")
    print("Structure information:")
    print(struc_copy)


def except_block(e):
    print(f"{str(e)}.")
    print("Please update pymatgen package to latest version better.")
    exit(1)


def hcp_gen(a, c, x_scale, y_scale, z_scale, species):
    lattice = Lattice.hexagonal(a, c)
    struc = Structure(
        lattice, [species for _ in range(2)], [[0.0, 0.0, 0.0], [2 / 3, 1 / 3, 1 / 2]]
    )
    struc_type = "HCP"
    supercell_gen(struc, x_scale, y_scale, z_scale, struc_type=struc_type)


def hcp_orth_gen(a, c, x_scale, y_scale, z_scale, species):
    lattice = Lattice.orthorhombic(a, a * math.sqrt(3), c)
    struc = Structure(
        lattice,
        [species for _ in range(4)],
        [
            [0.0, 0.0, 0.0],
            [1 / 2, 1 / 2, 0.0],
            [1 / 2, 1 / 3, 1 / 2],
            [0.0, 2 / 3, 1 / 2],
        ],
    )
    struc_type = "HCP Orthorhombic"
    supercell_gen(struc, x_scale, y_scale, z_scale, struc_type=struc_type)


def bcc_gen(a, x_scale, y_scale, z_scale, species):
    lattice = Lattice.cubic(a)
    struc = Structure(
        lattice, [species for _ in range(2)], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    )
    struc_tye = "BCC"
    supercell_gen(struc, x_scale, y_scale, z_scale, struc_type=struc_tye)


def b2_gen(a, x_scale, y_scale, z_scale, species_first, species_second):
    lattice = Lattice.cubic(a)
    struc = Structure(
        lattice,
        [f"{species_first}", f"{species_second}"],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )
    # struc = Structure.from_prototype(
    #     prototype="CsCl", species=[f"{species_first}", f"{species_second}"], a=a
    # )
    struc_type = "B2"
    supercell_gen(struc, x_scale, y_scale, z_scale, struc_type=struc_type)


def fcc_gen(a, x_scale, y_scale, z_scale, species):
    lattice = Lattice.cubic(a)
    struc = Structure(
        lattice,
        [species for _ in range(4)],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]],
    )
    struc_type = "FCC"
    supercell_gen(struc, x_scale, y_scale, z_scale, struc_type=struc_type)


def diamond_gen(a, x_scale, y_scale, z_scale, species):
    try:
        struc = Structure.from_prototype(
            prototype="diamond", species=[f"{species}"], a=a
        )
    except Exception as e:
        except_block(e)
    struc_type = "diamond"
    supercell_gen(struc, x_scale, y_scale, z_scale, struc_type=struc_type)


def rocksalt_gen(a, x_scale, y_scale, z_scale, species_first, species_second):
    try:
        struc = Structure.from_prototype(
            prototype="rocksalt", species=[f"{species_first}", f"{species_second}"], a=a
        )
    except Exception as e:
        except_block(e)
    struc_type = "rocksalt"
    supercell_gen(struc, x_scale, y_scale, z_scale, struc_type=struc_type)


def perovskite_gen(
    a, x_scale, y_scale, z_scale, species_first, species_second, species_third
):
    try:
        struc = Structure.from_prototype(
            prototype="perovskite",
            species=[f"{species_first}", f"{species_second}", f"{species_third}"],
            a=a,
        )
    except Exception as e:
        except_block(e)
    struc_type = "perovskite"
    supercell_gen(struc, x_scale, y_scale, z_scale, struc_type=struc_type)


def main():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(help="commands")

    hcp_args_list = [
        {"name": "a", "type": float, "help": "lattice constant a"},
        {"name": "c", "type": float, "help": "lattice constant c"},
        {"name": "x_scale", "type": int, "help": "x scale"},
        {"name": "y_scale", "type": int, "help": "y scale"},
        {"name": "z_scale", "type": int, "help": "z scale"},
        {"name": "species", "type": str, "help": "species"},
    ]
    cubic_args_list = hcp_args_list.copy()
    cubic_args_list.pop(1)

    b2_args_list = hcp_args_list.copy()
    b2_args_list.pop(1)
    b2_args_list.pop(-1)
    b2_args_list.append({"name": "species_first", "type": str, "help": "first species"})
    b2_args_list.append(
        {"name": "species_second", "type": str, "help": "second species"}
    )

    perovskite_args_list = b2_args_list.copy()
    perovskite_args_list.append(
        {"name": "species_third", "type": str, "help": "third species"}
    )

    commands = {
        "hcp": (hcp_gen, hcp_args_list),
        "hcp-o": (hcp_orth_gen, hcp_args_list),
        "bcc": (bcc_gen, cubic_args_list),
        "fcc": (fcc_gen, cubic_args_list),
        "b2": (b2_gen, b2_args_list),
        "diamond": (diamond_gen, cubic_args_list),
        "rocksalt": (rocksalt_gen, b2_args_list),
        "perovskite": (perovskite_gen, perovskite_args_list),
    }

    for command, (func, args_list) in commands.items():
        command_parser = subparsers.add_parser(
            name=command, help=f"{command} generation."
        )
        for arg_info in args_list:
            command_parser.add_argument(
                arg_info["name"], type=arg_info["type"], help=arg_info["help"]
            )
        command_parser.set_defaults(func=func)

    args = parser.parse_args()
    # 执行函数功能
    if hasattr(args, "func"):  # 检查是否存在处理函数
        if args.func == hcp_gen:
            hcp_gen(
                args.a, args.c, args.x_scale, args.y_scale, args.z_scale, args.species
            )
        elif args.func == hcp_orth_gen:
            hcp_orth_gen(
                args.a, args.c, args.x_scale, args.y_scale, args.z_scale, args.species
            )
        elif args.func == bcc_gen:
            bcc_gen(args.a, args.x_scale, args.y_scale, args.z_scale, args.species)
        elif args.func == b2_gen:
            b2_gen(
                args.a,
                args.x_scale,
                args.y_scale,
                args.z_scale,
                args.species_first,
                args.species_second,
            )
        elif args.func == fcc_gen:
            fcc_gen(args.a, args.x_scale, args.y_scale, args.z_scale, args.species)
        elif args.func == diamond_gen:
            diamond_gen(args.a, args.x_scale, args.y_scale, args.z_scale, args.species)
        elif args.func == rocksalt_gen:
            rocksalt_gen(
                args.a,
                args.x_scale,
                args.y_scale,
                args.z_scale,
                args.species_first,
                args.species_second,
            )
        elif args.func == perovskite_gen:
            perovskite_gen(
                args.a,
                args.x_scale,
                args.y_scale,
                args.z_scale,
                args.species_first,
                args.species_second,
                args.species_third,
            )


if __name__ == "__main__":
    main()
