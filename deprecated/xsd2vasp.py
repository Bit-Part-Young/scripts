"""构型文件格式转换: Material Studio XSD 与 VASP"""

import argparse

from ase.io import read, write


def xsd2vasp(xsd_file: str, output_file: str = None):
    atoms = read(xsd_file, format="xsd")

    if not output_file:
        output_file = f"{xsd_file.rsplit('.', 1)[0]}.vasp"

    write(output_file, images=atoms, format="vasp", direct=True, sort=True)


def vasp2xsd(vasp_file: str, output_file: str = None):
    atoms = read(vasp_file, format="vasp")

    if not output_file:
        if vasp_file == "POSCAR":
            fn = atoms.get_chemical_formula()
        elif vasp_file.endswith((".vasp", ".poscar")):
            fn = vasp_file.rsplit(".", 1)[0]

        output_file = f"{fn}.xsd"

    write(output_file, images=atoms, format="xsd")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Configuration file format conversion between Material Studio XSD and VASP."
    )
    parser.add_argument(
        "--option",
        default="xsd2vasp",
        choices=["xsd2vasp", "vasp2xsd"],
        required=True,
        help="Conversion option.",
    )
    parser.add_argument("-i", "--input", required=True, help="Input file.")
    parser.add_argument("-o", "--output", help="Output file, optional.")

    args = parser.parse_args()

    if args.option == "xsd2vasp":
        xsd2vasp(args.input, args.output)
    elif args.option == "vasp2xsd":
        vasp2xsd(args.input, args.output)
