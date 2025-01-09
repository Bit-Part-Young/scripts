"""FCC (100), (110), (111) 表面模型构建"""

import os

from ase.build import fcc100, fcc110, fcc111
from ase.io import read, write


def get_fcc_surface(
    symbol: str,
    a: float,
    surface_type: str,
    slab_size: tuple[int] = (1, 1, 6),
    vacuum: float = 5.0,
    output_fn: str = "-",
    output_format: str = "vasp",
):
    if surface_type == "100":
        atoms = fcc100(
            symbol=symbol,
            a=a,
            size=slab_size,
            vacuum=vacuum,
            orthogonal=True,
        )
    elif surface_type == "110":
        atoms = fcc110(
            symbol=symbol,
            a=a,
            size=slab_size,
            vacuum=vacuum,
            orthogonal=True,
        )
    elif surface_type == "111":
        atoms = fcc111(
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

    symbol = "Cu"
    a = 3.615

    file_100 = f"{symbol}_100.lmp"
    file_110 = f"{symbol}_110.lmp"
    file_111 = f"{symbol}_111.lmp"

    # (100) 表面
    get_fcc_surface(
        symbol=symbol,
        a=a,
        surface_type="100",
        slab_size=(1, 1, 6),
        output_fn=file_100,
    )

    # (110) 表面
    get_fcc_surface(
        symbol=symbol,
        a=a,
        surface_type="110",
        slab_size=(1, 1, 6),
        output_fn=file_110,
    )

    # (111) 表面
    get_fcc_surface(
        symbol=symbol,
        a=a,
        surface_type="111",
        slab_size=(1, 2, 6),
        output_fn=file_111,
    )

    surface_folder = "surfaces-ase"
    os.makedirs(surface_folder, exist_ok=True)
    os.system(f"mv {file_100} {file_110} {file_111} {surface_folder}")

    """
    atoms = read(
        file_110,
        format="lammps-data",
        units="metal",
        atom_style="atomic",
        Z_of_type={
            1: atomic_numbers["Cu"],
        },
    )
    """
