#!/usr/bin/env python3

"""生成含位移的超胞，用于 phonopy 计算"""

import argparse

import phonopy
from phonopy.interface.calculator import (
    write_crystal_structure,
    write_supercells_with_displacements,
)
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core.structure import Structure
from pymatgen.io.phonopy import get_phonopy_structure

structure_fn_dict = {
    "vasp": ["PPOSCAR", "BPOSCAR", "SPOSCAR"],
    "lammps": ["Punitcell", "Bunitcell", "supercell"],
}


def unit2primitive(
    structure_fn: str = "POSCAR", interface_mode: str = "vasp"
) -> PhonopyAtoms:
    structure = Structure.from_file(structure_fn)

    primitive_structure = structure.to_primitive()

    write_crystal_structure(
        filename=structure_fn_dict[interface_mode][1],
        cell=get_phonopy_structure(structure),
        interface_mode=interface_mode,
    )

    return get_phonopy_structure(primitive_structure)


def disp_supercell_generation(
    structure_fn: str = "POSCAR",
    interface_mode: str = "vasp",
    cellsize: int = 4,
):
    """生成含位移的超胞"""

    unitcell = unit2primitive(structure_fn, interface_mode)

    phonon = phonopy.load(
        unitcell=unitcell,
        primitive_matrix="auto",
        supercell_matrix=[cellsize, cellsize, cellsize],
        calculator=interface_mode,
    )

    phonon.generate_displacements()
    phonon.save("phonopy_disp.yaml")

    write_crystal_structure(
        filename=structure_fn_dict[interface_mode][0],
        cell=phonon.primitive,
        interface_mode=interface_mode,
    )

    write_crystal_structure(
        filename=structure_fn_dict[interface_mode][2],
        cell=phonon.supercell,
        interface_mode=interface_mode,
    )

    write_supercells_with_displacements(
        supercell=phonon.supercell,
        cells_with_disps=phonon.supercells_with_displacements,
        interface_mode=interface_mode,
    )

    if interface_mode == "vasp":
        print(f"\nPPOSCAR, BPOSCAR, SPOSCAR and POSCAR-XXX generated.")
    elif interface_mode == "lammps":
        print(f"\nPunitcell, Bunitcell, supercell and supercell-XXX generated.")
    print(f"phonopy_disp.yaml generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate supercells with displacements for phonopy calculation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )

    parser.add_argument(
        "structure_fn",
        type=str,
        default="POSCAR",
        help="structure filename",
    )

    parser.add_argument(
        "--cellsize",
        type=int,
        default=4,
        const=4,
        nargs="?",
        metavar="N",
        help="cellsize",
    )

    parser.add_argument(
        "--interface",
        type=str,
        default="vasp",
        const="vasp",
        nargs="?",
        choices=["vasp", "lammps"],
        help="interface mode.",
    )

    args = parser.parse_args()

    disp_supercell_generation(
        structure_fn=args.structure_fn,
        interface_mode=args.interface,
        cellsize=args.cellsize,
    )
