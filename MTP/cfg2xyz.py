"""
Purpose:
    Convert mtp input file format to xyz.
Ref:
    dpdata: https://github.com/deepmodeling/dpdata
Run:
    python mtp2xyz.py train.cfg Symbol1 Symbol2 Symbol3 ...
"""

import os
import sys
from collections import defaultdict
from typing import Dict

import numpy as np
from ase.atoms import Atoms


def load_cfg(cfg_fn: str, type_to_symbol: dict[int, str]) -> list[Atoms]:
    frames = []
    with open(cfg_fn) as f:
        line = "chongchongchong!"
        while line:
            line = f.readline()
            if "BEGIN_CFG" in line:
                cell = np.zeros((3, 3))
            if "Size" in line:
                line = f.readline()
                natoms = int(line.split()[0])
                positions = np.zeros((natoms, 3))
                forces = np.zeros((natoms, 3))
                energies = np.zeros(natoms)
                symbols = ["X"] * natoms
            if "Supercell" in line:
                for fields_key in range(3):
                    line = f.readline()
                    for j in range(3):
                        cell[fields_key, j] = float(line.split()[j])
            if "AtomData" in line:
                fields_kv_dict = defaultdict(int)
                for fields_key, fields_value in enumerate(line.split()[1:]):
                    fields_kv_dict[fields_value] = fields_key

                for _ in range(natoms):
                    line = f.readline()
                    fields_content = line.split()
                    atom_index = int(fields_content[fields_kv_dict["id"]]) - 1
                    symbols[atom_index] = type_to_symbol[
                        int(fields_content[fields_kv_dict["type"]])
                    ]
                    positions[atom_index] = [
                        float(fields_content[fields_kv_dict[key]])
                        for key in ["cartes_x", "cartes_y", "cartes_z"]
                    ]
                    forces[atom_index] = [
                        float(fields_content[fields_kv_dict[key]])
                        for key in ["fx", "fy", "fz"]
                    ]

                atoms = Atoms(symbols=symbols, cell=cell, positions=positions, pbc=True)
                if fields_kv_dict["fx"] != 0:
                    atoms.info["forces"] = forces

            if "Energy" in line and "Weight" not in line:
                line = f.readline()
                atoms.info["energy"] = float(line.split()[0])
            if "PlusStress" in line:
                line = f.readline()
                # MTP cfg stress 分量顺序 xx yy zz yz xz xy
                plusstress = np.array(list(map(float, line.split())))
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
                frames.append(atoms)

    return frames


def dump_xyz(frames: list[Atoms]):

    xyz_string = ""
    nframes = len(frames)
    for atoms in frames:
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
        forces = atoms.info["forces"]
        for i in range(len(atoms)):
            xyz_string += "{:2} {:>15.8f} {:>15.8f} {:>15.8f} {:>15.8f} {:>15.8f} {:>15.8f}\n".format(
                symbols[i], *positions[i], *forces[i]
            )

    os.makedirs("XYZ", exist_ok=True)
    with open(os.path.join("XYZ", "mtp2xyz.xyz"), "w") as f:
        f.write(xyz_string)


if __name__ == "__main__":
    type_to_symbol = {i: s for i, s in enumerate(sys.argv[2:])}
    frames = load_cfg(sys.argv[1], type_to_symbol)
    dump_xyz(frames)
