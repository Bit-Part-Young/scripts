#!/usr/bin/env python3

"""
将 NEP xyz 文件转换成 MTP cfg

reference: https://github.com/hu-yanxiao/SUS2-MLIP/blob/main/python_tools/xyz2cfg.py
"""

import argparse

from ase.io import read


def xyz2cfg(
    xyz_fn: str,
    symbols_map: dict[str, int],
    cfg_fn: str = "output.cfg",
):
    """将 NEP xyz 文件转换为 MTP cfg"""

    atoms_list = read(xyz_fn, index=":", format="extxyz")

    ff = open(cfg_fn, "w")
    for atoms in atoms_list:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        plusstress = atoms.info["virial"]

        chemical_symbols = atoms.get_chemical_symbols()
        natoms = len(chemical_symbols)

        cell = atoms.get_cell()
        positioins = atoms.get_positions()

        ff.write("""BEGIN_CFG\n""")
        ff.write(""" Size\n""")
        ff.write("""  {:6}\n""".format(natoms))
        ff.write(""" Supercell \n""")
        ff.write(
            """{:12.8f} {:12.8f} {:12.8f}\n""".format(
                cell[0, 0], cell[0, 1], cell[0, 2]
            )
        )
        ff.write(
            """{:12.8f} {:12.8f} {:12.8f}\n""".format(
                cell[1, 0], cell[1, 1], cell[1, 2]
            )
        )
        ff.write(
            """{:12.8f} {:12.8f} {:12.8f}\n""".format(
                cell[2, 0], cell[2, 1], cell[2, 2]
            )
        )
        ff.write(
            """ AtomData:  id   type         cartes_x         cartes_y         cartes_z         fx         fy         fz\n"""
        )
        for atom_index in range(natoms):
            ff.write(
                """        {:6} {:6} {:12.8f}   {:12.8f}   {:12.8f}   {:12.8f}   {:12.8f}   {:12.8f}\n""".format(
                    atom_index + 1,
                    symbols_map[chemical_symbols[atom_index]],
                    positioins[atom_index, 0],
                    positioins[atom_index, 1],
                    positioins[atom_index, 2],
                    forces[atom_index, 0],
                    forces[atom_index, 1],
                    forces[atom_index, 2],
                )
            )
        ff.write(""" Energy \n""")
        ff.write(f"""     {energy:12.8f} \n""")
        ff.write(
            """ PlusStress:  xx          yy          zz          yz          xz          xy\n"""
        )
        ff.write(
            "     {:12.8f} {:12.8f} {:12.8f} {:12.8f} {:12.8f} {:12.8f}\n".format(
                plusstress[0, 0],
                plusstress[1, 1],
                plusstress[2, 2],
                plusstress[1, 2],
                plusstress[0, 2],
                plusstress[0, 1],
            )
        )
        ff.write("""END_CFG\n""")
    ff.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert NEP xyz to cfg.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("xyz_fn", type=str, default="train.xyz", help="xyz filename")
    parser.add_argument("cfg_fn", type=str, default="output.cfg", help="cfg filename")
    parser.add_argument(
        "-ess",
        type=str,
        nargs="+",
        required=True,
        help="element symbol sequences, e.g. Ti Al Nb",
    )

    args = parser.parse_args()
    xyz_fn = args.xyz_fn
    cfg_fn = args.cfg_fn
    ess = args.ess

    symbols_map = {symbol: index for index, symbol in enumerate(ess)}

    xyz2cfg(xyz_fn=xyz_fn, symbols_map=symbols_map, cfg_fn=cfg_fn)

    print(f"Convert {xyz_fn} to {cfg_fn}.")
