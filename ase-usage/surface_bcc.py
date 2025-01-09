"""FCC (100), (110), (111) 表面模型构建"""

import os

from ase.build import bcc100, bcc110, bcc111
from ase.io import read, write


def get_bcc_surface(
    symbol: str,
    a: float,
    surface_type: str,
    slab_size: tuple[int] = (1, 1, 6),
    vacuum: float = 5.0,
    output_fn: str = "-",
    output_format: str = "vasp",
):
    if surface_type == "100":
        atoms = bcc100(
            symbol=symbol,
            a=a,
            size=slab_size,
            vacuum=vacuum,
            orthogonal=True,
        )
    elif surface_type == "110":
        atoms = bcc110(
            symbol=symbol,
            a=a,
            size=slab_size,
            vacuum=vacuum,
            orthogonal=True,
        )
    elif surface_type == "111":
        atoms = bcc111(
            symbol=symbol,
            a=a,
            size=slab_size,
            vacuum=vacuum,
            orthogonal=True,
        )

    if output_format == "vasp":
        write(
            output_fn,
            atoms,
            format="vasp",
            direct=True,
            sort=True,
        )
    elif output_format == "lammps-data":
        write(
            output_fn,
            atoms,
            format="lammps-data",
            units="metal",
            atom_style="atomic",
            specorder=[symbol],
        )


if __name__ == "__main__":

    symbol = "Nb"
    a = 3.32

    file_100 = f"{symbol}_100.vasp"
    file_110 = f"{symbol}_110.vasp"
    file_111 = f"{symbol}_111.vasp"

    # (100) 表面
    get_bcc_surface(
        symbol=symbol,
        a=a,
        surface_type="100",
        slab_size=(1, 1, 6),
        vacuum=10.0,
        output_fn=file_100,
    )

    # (110) 表面
    get_bcc_surface(
        symbol=symbol,
        a=a,
        surface_type="110",
        slab_size=(1, 2, 6),
        vacuum=10.0,
        output_fn=file_110,
    )

    # (111) 表面
    get_bcc_surface(
        symbol=symbol,
        a=a,
        surface_type="111",
        slab_size=(1, 2, 6),
        vacuum=10.0,
        output_fn=file_111,
    )

    surface_folder = "surfaces-ase"
    os.makedirs(surface_folder, exist_ok=True)
    os.system(f"mv {file_100} {file_110} {file_111} {surface_folder}")
