#!/usr/bin/env python3

"""构型格式转换"""

import argparse

from ase.io import read, write


def posconv(
    input_file: str,
    output_file: str,
):
    """构型格式转换"""

    atoms = read(input_file)

    if (output_file == "POSCAR") or output_file.endswith("vasp"):

        write(
            output_file,
            images=atoms,
            format="vasp",
            direct=True,
            sort=True,
        )

    elif output_file.endswith("xsd"):
        write(
            output_file,
            images=atoms,
            format="xsd",
        )

    elif output_file.endswith("xyz"):
        write(
            output_file,
            images=atoms,
            format="xyz",
        )

    elif (input_file == "OUTCAR") & output_file.endswith("xyz"):
        write(
            output_file,
            images=atoms,
            format="extxyz",
        )

    print(f"\nConvert {input_file} to {output_file} done!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Structure file format convert. Support most ASE recognized formats.",
        epilog="Author: SLY.",
        allow_abbrev=True,
    )

    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        required=True,
        choices=[
            "POSCAR",
            "CONTCAR",
            "OUTCAR",
            "*.vasp",
            "*.xsd",
            "*.xyz",
        ],
        help="Input file.",
    )

    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        choices=[
            "POSCAR",
            "*.vasp",
            "*.xsd",
            "*.xyz",
        ],
        help="Output file.",
    )

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    posconv(
        input_file=input_file,
        output_file=output_file,
    )
