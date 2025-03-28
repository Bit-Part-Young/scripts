#!/usr/bin/env python3

"""将 MTP cfg 文件转换成 NEP xyz"""

import argparse
import os

import numpy as np
from ase.atoms import Atoms


def load_cfg(cfg_fn: str, symbols_map: dict[int, str]) -> list[Atoms]:
    """导入 cfg 文件，返回 list[Atoms]"""

    atoms_list = []
    with open(cfg_fn) as f:
        line = "chongchongchong!"
        while line:
            line = f.readline()
            if "Size" in line:
                line = f.readline()
                natoms = int(line.split()[0])

            if "Supercell" in line:
                cell = []
                for _ in range(3):
                    line = f.readline()
                    cell.append([float(x) for x in line.split()])
                cell = np.array(cell)

            if "AtomData" in line:
                atomdata_list = []
                for _ in range(natoms):
                    line = f.readline()
                    atomdata_list.append([float(x) for x in line.split()])

                atomdata_array = np.array(atomdata_list)

                symbols_int = atomdata_array[:, 1].tolist()
                symbols = [symbols_map[symbol_int] for symbol_int in symbols_int]

                positions = atomdata_array[:, 2:5]
                forces = atomdata_array[:, 5:8]

                atoms = Atoms(symbols=symbols, cell=cell, positions=positions, pbc=True)
                atoms.arrays["forces"] = forces

            if "Energy" in line and "Weight" not in line:
                line = f.readline()
                atoms.info["energy"] = float(line.split()[0])

            if "PlusStress" in line:
                line = f.readline()
                # MTP cfg stress 分量顺序 xx yy zz yz xz xy
                plusstress = list(map(float, line.split()))
                # MTP cfg 中的 virial 对应与 NEP xyz 中的 varial
                virial = [
                    plusstress[0],
                    plusstress[5],
                    plusstress[4],
                    plusstress[5],
                    plusstress[1],
                    plusstress[3],
                    plusstress[4],
                    plusstress[3],
                    plusstress[2],
                ]
                atoms.info["virial"] = virial

            if "END_CFG" in line:
                atoms_list.append(atoms)

    return atoms_list


def dump_xyz(atoms_list: list[Atoms], xyz_fn: str = "mtp2xyz.xyz"):
    """将 list[Atoms] 转换为 xyz 文件"""

    xyz_string = ""
    for atoms in atoms_list:
        xyz_string += str(len(atoms)) + "\n"
        xyz_string += "config_type=cfg2xyz "
        xyz_string += "energy=" + str(atoms.info["energy"]) + " "
        xyz_string += 'pbc="T T T" '
        xyz_string += 'virial="' + " ".join(list(map(str, atoms.info["virial"]))) + '" '
        xyz_string += (
            'Lattice="' + " ".join(list(map(str, atoms.get_cell().reshape(-1)))) + '" '
        )
        xyz_string += "Properties=species:S:1:pos:R:3:forces:R:3\n"

        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        forces = atoms.arrays["forces"]
        for i in range(len(atoms)):
            xyz_string += "{:2} {:>15.8f} {:>15.8f} {:>15.8f} {:>15.8f} {:>15.8f} {:>15.8f}\n".format(
                symbols[i], *positions[i], *forces[i]
            )

    with open(xyz_fn, "w") as f:
        f.write(xyz_string)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert MTP cfg to NEP xyz.")

    parser.add_argument(
        "cfg_fn",
        type=str,
        default="train.cfg",
        help="MTP cfg filename",
    )

    parser.add_argument(
        "xyz_fn",
        type=str,
        default="mtp2xyz.xyz",
        help="NEP xyz filename",
    )

    parser.add_argument(
        "-ess",
        type=str,
        nargs="+",
        required=True,
        help="element symbol sequences, e.g. Ti Al Nb",
    )

    args = parser.parse_args()
    cfg_fn = args.cfg_fn
    xyz_fn = args.xyz_fn
    ess = args.ess

    symbols_map = {i: s for i, s in enumerate(ess)}
    atoms_list = load_cfg(cfg_fn, symbols_map)
    dump_xyz(atoms_list, xyz_fn)

    print(f"Convert {cfg_fn} to {xyz_fn}.")
