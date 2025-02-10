#!/usr/bin/env python3

"""根据计算类型生成 VASP INCAR 文件"""

import argparse


def incar_template(
    calculation_type: str,
    system: str,
    icharg: int,
    encut: int,
    ismear: int,
    sigma: float,
    ediff: float,
    nelm: int,
    nsw: int,
    ibrion: int,
    isif: int,
    ediffg: float,
):
    """INCAR 文件模板"""

    with open("INCAR", "w") as f:
        f.write("Global Parameters\n")
        f.write(f"SYSTEM  = {system}\n")
        f.write("ISTART  = 0\n")
        f.write(f"ICHARG  = {icharg}\n")
        f.write("ISPIN   = 1\n")
        f.write("LCHARG  = .FALSE.\n")
        f.write("LWAVE   = .FALSE.\n")

        if calculation_type == "phonon":
            f.write("ADDGRID = .TRUE.\n")

        f.write("PREC    = Accurate\n")
        f.write("LREAL   = .FALSE.\n")
        f.write(f"ENCUT   = {encut}\n\n")

        f.write("Electronic Relaxation\n")
        f.write("ALGO    = Fast\n")
        f.write(f"ISMEAR  = {ismear}\n")
        f.write(f"SIGMA   = {sigma}\n")
        f.write(f"EDIFF   = {ediff:.0E}\n")
        f.write(f"NELM    = {nelm}\n")
        f.write("NELMIN  = 6\n\n")

        f.write("Ionic Relxation\n")
        f.write(f"NSW     = {nsw}\n")
        f.write(f"IBRION  = {ibrion}\n")
        f.write(f"ISIF    = {isif}\n")

        if calculation_type in ["relax", "phonon"]:
            f.write(f"EDIFFG  = {ediffg:.0E}\n")

        elif calculation_type == "dos":
            f.write("\nLORBIT  = 11\n")

    print(f"INCAR {calculation_type} template generated.\n")


def incar_generation(
    calculation_type: str = "scf",
    system: str = "Static",
    icharg: int = 2,
    encut: int = 400,
    ismear: int = -5,
    sigma: float = 0.05,
    ediff: float = 1e-6,
    nelm: int = 90,
    nsw: int = 0,
    ibrion: int = -1,
    isif: int = 2,
    ediffg: float = -1e-2,
):
    """生成 INCAR 文件"""

    incar_dict = {
        "system": system,
        "icharg": icharg,
        "encut": encut,
        "ismear": ismear,
        "sigma": sigma,
        "ediff": ediff,
        "nelm": nelm,
        "nsw": nsw,
        "ibrion": ibrion,
        "isif": isif,
        "ediffg": ediffg,
    }

    if calculation_type == "relax":
        incar_dict.update(
            {
                "system": "Relaxation",
                "nsw": 100,
                "ismear": 1,
                "ibrion": 2,
            }
        )
    elif calculation_type == "dos":
        incar_dict.update(
            {
                "system": "DOS",
                "icharg": 11,
            }
        )
    elif calculation_type == "band":
        incar_dict.update(
            {
                "system": "Band",
                "icharg": 11,
                "ismear": 0,
            }
        )
    elif calculation_type == "phonon":
        incar_dict.update(
            {
                "system": "Phonon",
                "ismear": 0,
                "ediff": 1e-8,
                "nsw": 1,
                "ibrion": 8,
            }
        )

    incar_template(calculation_type=calculation_type, **incar_dict)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate VASP input files.",
        epilog="Author: SLY.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=True,
    )

    parser.add_argument(
        "-ct",
        "--calculation_type",
        nargs="?",
        const="scf",
        default="scf",
        type=str,
        help="Calculation type",
    )

    group = parser.add_argument_group("INCAR tags")
    group.add_argument("--encut", nargs="?", const=400, default=400, type=int, help="ENCUT tag")
    group.add_argument("--ediff", nargs="?", const=1e-6, default=1e-6, type=float, help="EDIFF tag")
    group.add_argument("--nelm", nargs="?", const=90, default=90, type=int, help="NELM tag")
    group.add_argument("--isif", nargs="?", const=2, default=2, type=int, help="ISIF tag")

    args = parser.parse_args()

    calculation_type = args.calculation_type
    encut = args.encut
    ediff = args.ediff
    nelm = args.nelm
    isif = args.isif

    incar_generation(
        calculation_type=calculation_type,
        encut=encut,
        ediff=ediff,
        nelm=nelm,
        isif=isif,
    )
