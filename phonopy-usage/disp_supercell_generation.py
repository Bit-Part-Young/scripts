#!/usr/bin/env python3

"""
生成 phonopy 有限位移超胞（VASP 或 LAMMPS 格式），用于 DFT 计算或用 LAMMPS 校验势函数的声子谱

原胞由 phonopy 自动识别（primitive_matrix = auto，spglib 标准设定，hcp 为 120°），与 band_*.conf 的 k 路径一致。
输出：phonopy_disp.yaml；vasp：PPOSCAR、SPOSCAR、POSCAR-XXX；lammps：Punitcell、supercell、supercell-XXX

LAMMPS 流程：对每个 supercell-XXX 做单点计算，按 phonopy 文档的格式 dump 力（dump_modify sort id），然后
    phonopy-init --lammps -f lammps_forces_*.0
    phonopy --config band_hcp.conf -s
"""

import argparse

import phonopy
from ase.io import read
from phonopy.interface.calculator import write_crystal_structure, write_supercells_with_displacements
from phonopy.structure.atoms import PhonopyAtoms

OUTPUT_NAMES = {"vasp": ("PPOSCAR", "SPOSCAR"), "lammps": ("Punitcell", "supercell")}


def disp_supercell_generation(
    structure_fn: str = "POSCAR",
    interface_mode: str = "vasp",
    dim: tuple[int, int, int] = (4, 4, 4),
    distance: float = 0.01,
):
    atoms = read(structure_fn)
    unitcell = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.cell[:],
        scaled_positions=atoms.get_scaled_positions(),
    )

    phonon = phonopy.Phonopy(
        unitcell, supercell_matrix=dim, primitive_matrix="auto", calculator=interface_mode
    )
    phonon.generate_displacements(distance=distance)
    phonon.save("phonopy_disp.yaml")

    primitive_fn, supercell_fn = OUTPUT_NAMES[interface_mode]
    write_crystal_structure(primitive_fn, phonon.primitive, interface_mode=interface_mode)
    write_crystal_structure(supercell_fn, phonon.supercell, interface_mode=interface_mode)
    write_supercells_with_displacements(
        interface_mode, phonon.supercell, phonon.supercells_with_displacements
    )

    n = len(phonon.supercells_with_displacements)
    print(f"phonopy_disp.yaml, {primitive_fn}, {supercell_fn} and {n} displaced supercells generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate supercells with displacements for phonopy calculation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Author: SLY.",
    )
    parser.add_argument("structure_fn", nargs="?", default="POSCAR", help="structure file (any ASE format)")
    parser.add_argument("-d", "--dim", nargs=3, type=int, default=[4, 4, 4], metavar="N", help="supercell size")
    parser.add_argument("--distance", type=float, default=0.01, help="displacement distance (Å)")
    parser.add_argument("--interface", default="vasp", choices=list(OUTPUT_NAMES), help="interface mode")
    args = parser.parse_args()

    disp_supercell_generation(args.structure_fn, args.interface, args.dim, args.distance)
