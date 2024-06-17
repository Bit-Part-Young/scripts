import argparse

import numpy as np
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp.inputs import Poscar


def pos2cif(poscar_filename="POSCAR"):
    """
    Transform poscar file to cif file.
    little complicated.

    To do:
        - get the right cif filename prefix. Done
        - get the right symmetry of the structure and save primitive cif.
    """
    pos = Poscar.from_file(poscar_filename)
    struc = pos.structure
    species = pos.site_symbols
    natoms = pos.natoms

    # get reduced atoms number list
    gcd_result = np.gcd.reduce(natoms)
    natoms_reduced = ((np.array(natoms) / gcd_result).astype(int)).tolist()
    natoms_reduced_str = [str(i) for i in natoms_reduced]

    reduced_formula = "".join(
        [i + j if j != "1" else i for i, j in zip(species, natoms_reduced_str)]
    )
    # print(reduced_formula)

    cif_fn = f"{reduced_formula}.cif"
    CifWriter(struc).write_file(cif_fn)
    # CifWriter(struc, symprec=0.01).write_file(cif_fn)

    print(f"{poscar_filename} is transformed to {cif_fn}.")


if __name__ == "__main__":
    args = argparse.ArgumentParser("transform poscar to cif file.")
    args.add_argument("-p", "--poscar", type=str, help="poscar file name")

    args = args.parse_args()

    pos2cif(poscar_filename=args.poscar)
