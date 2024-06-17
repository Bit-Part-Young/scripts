import argparse

from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Poscar


def cif2pos(cif_filename, primitive=False):
    """
    Transform cif file to poscar file.
    """
    cp = CifParser(cif_filename)
    struct = cp.parse_structures(primitive=primitive)[0]

    poscar_filename = f"{cif_filename[:-4]}.vasp"
    Poscar(struct).write_file(poscar_filename)

    if not primitive:
        print(f"{cif_filename} is transformed to {poscar_filename} with unit cell.")
    else:
        print(
            f"{cif_filename} is transformed to {poscar_filename} with primitive cell."
        )


if __name__ == "__main__":
    args = argparse.ArgumentParser("transform cif to poscar file.")
    args.add_argument("-c", "--cif", type=str, help="cif filename")
    args.add_argument(
        "-p", "--primitive", type=int, default=0, help="whether primitive"
    )

    args = args.parse_args()

    cif2pos(cif_filename=args.cif, primitive=args.primitive)
