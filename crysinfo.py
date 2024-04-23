# reference: https://github.com/Shibu778/structure2symmetry/blob/main/structure2symmetry/symmetry.py

import os

import click
import numpy as np
from ase.io import read
from spglib import get_spacegroup, get_symmetry_dataset


def has_inversion(n_rotations):
    """Returns true if the n_rotations list has 3x3 operation for inversion.

    Parameters
    ----------
    n_rotations : list of 3x3 matrices
        List of n rotation matrix (Point Group Symmetry)

    Returns
    -------
    bool
        True for if n_rotations has inversion operator.
    """

    inversion = -np.identity(3, dtype=int)
    return np.any([np.all(rotation == inversion) for rotation in n_rotations])


def write_dataset(filename, dataset):
    """Takes dataset and write it in a file beautifully.

    Check the following reference for details about various terms.
    Ref - https://spglib.readthedocs.io/en/latest/python-spglib.html

    Parameters
    ----------
    filename : str or os.path
        File path with filename where the symmetry dataset should be written
    dataset : dict
        Contains all the information to be written into the output file

    """

    lines = []
    lines.append(
        "Check the following reference for details about various terms. \nRef - https://spglib.readthedocs.io/en/latest/python-spglib.html"
    )
    lines.append("\n\nFile Read : " + str(dataset["filename_read"]) + "\n")
    lines.append("\n\nSpace Group Number : " + str(dataset["number"]) + "\n")
    lines.append("Point Group : " + str(dataset["pointgroup"]) + "\n")
    lines.append(
        "Space Group (International) : " + str(dataset["international"]) + "\n"
    )
    lines.append("Space Group (Schoenflies) : " + str(dataset["schoenflies"]) + "\n")
    lines.append("Has inversion symmetry : " + str(dataset["has_inversion"]) + "\n")
    lines.append("Symmetry Precision : " + str(dataset["tolerance"]) + "\n")
    lines.append("Angle Tolerance : " + str(dataset["angle_tolerance"]) + "\n")
    lines.append("Hall Number : " + str(dataset["hall_number"]) + "\n")
    lines.append("Hall Symbol : " + str(dataset["hall"]) + "\n")
    lines.append("Choice : " + str(dataset["choice"]) + "\n")
    lines.append("Transformation Matrix : \n" + str(dataset["transformation_matrix"]))
    lines.append("\n\nOrigin Shift : " + str(dataset["origin_shift"]))
    lines.append("\n\nRotations : \n" + str(dataset["rotations"]))
    lines.append("\n\nTranslations : \n" + str(dataset["translations"]))
    lines.append("\n\nWickoff Letters : \n" + str(dataset["wyckoffs"]))
    lines.append(
        "\n\nSite Symmetry Symbols : \n" + str(dataset["site_symmetry_symbols"])
    )
    lines.append(
        "\n\nCrystallographic Orbits : \n" + str(dataset["crystallographic_orbits"])
    )
    lines.append("\n\nEquivalent Atoms : \n" + str(dataset["equivalent_atoms"]))
    lines.append("\n\nPrimitive Lattice : \n" + str(dataset["primitive_lattice"]))
    lines.append("\n\nMapping to Primitive : \n" + str(dataset["mapping_to_primitive"]))
    lines.append("\n\nStandard Lattice : \n" + str(dataset["std_lattice"]))
    lines.append("\n\nStandard Positions : \n" + str(dataset["std_positions"]))
    lines.append(
        "\n\nStandard Rotation Matrix : \n" + str(dataset["std_rotation_matrix"])
    )
    lines.append(
        "\n\nStandard Mapping to Primitive : \n"
        + str(dataset["std_mapping_to_primitive"])
    )

    with open(filename, "w") as f:
        f.writelines(lines)


@click.command()
@click.option(
    "--filepath",
    default="./POSCAR",
    help="File path absolute or relative including the filename",
)
@click.option("--tolerance", default=1e-5, help="spglib symmetry tolerance")
@click.option("--angle_tolerance", default=0.01, help="spglib angle tolerance")
@click.option(
    "--destination",
    default="./",
    help="Destination folder where the symmetry information will be written",
)
@click.option(
    "--savefile",
    default="symmetry.info",
    help="Name of the file where symmetry information will be written",
)
def struct2symm(
    filepath,
    tolerance,
    angle_tolerance,
    destination="./",
    savefile="symmetry.info",
):
    """Returns a dictionary containing symmetry information.

    It uses spglib to determine the symmetry and ase to read the structure file.

    Supported structure files same as mentioned in https://wiki.fysik.dtu.dk/ase/ase/io/io.html .

    Parameters
    ----------
    filepath : str
        Absolute or relative filepath of the structure file
    tolerance : float
        Spglib symmetry precission
    angle_tolerance : float
        spglib angle precission
    destination : str or os.path
        Destination folder where the symmetry information will be written
    savefile : str
        Name of the file where symmetry information will be written

    Returns
    -------
    dataset : dictionary
        It contains all the information returned by dataset from spglib.get_symmetry_dataset and
        another key with has_inversion telling whether the inversion center is present or not.
    """

    struct = read(filepath)

    lattice = struct.get_cell().array
    positions = struct.get_scaled_positions(wrap=False)
    numbers = struct.get_atomic_numbers()

    spgcell = (lattice, positions, numbers)

    dataset = get_symmetry_dataset(
        spgcell, symprec=tolerance, angle_tolerance=angle_tolerance
    )

    dataset["has_inversion"] = has_inversion(dataset["rotations"])
    spacegroup_schoen = get_spacegroup(
        spgcell, symprec=tolerance, angle_tolerance=angle_tolerance, symbol_type=1
    )
    dataset["schoenflies"] = spacegroup_schoen
    dataset["tolerance"] = tolerance
    dataset["angle_tolerance"] = angle_tolerance
    dataset["filename_read"] = filepath

    filename = os.path.join(destination, savefile)
    write_dataset(filename, dataset)

    print("Point Group Symmetry : ", dataset["pointgroup"])
    print("Space Group Symmetry (International)", dataset["international"])
    print("Space Group Symmetry (Schoenflies)", dataset["schoenflies"])

    return dataset


if __name__ == "__main__":
    struct2symm()
#     # filepath = "../tests/POSCAR"
#     filepath = "../tests/SiC.cif"
#     destination = "../tests/"
#     tolerance = 1e-5
#     angle_tolerance = 0.01
#     dataset = struct2symm(filepath, tolerance, angle_tolerance, destination=destination)
#     # print(dataset)
