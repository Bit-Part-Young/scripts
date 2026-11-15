"""FCC (100), (110), (111) 表面模型构建（ase.build.fcc100/110/111）

注意：ASE 的 slab_size[2] 为原子层数，vacuum 为单侧真空厚度（总真空 = 2 * vacuum）；
取向与 atomsk / structure-scripts/build_surface.py 不同，需与 atomsk 一致时请用 build_surface.py。
"""

import shutil
from pathlib import Path

from ase.build import fcc100, fcc110, fcc111
from ase.io import write


def get_fcc_surface(
    symbol: str,
    a: float,
    surface_type: str,
    slab_size: tuple[int, int, int] = (1, 1, 6),
    vacuum: float = 5.0,
    output_fn: str = "-",
    output_format: str = "vasp",
):
    builders = {"100": fcc100, "110": fcc110, "111": fcc111}
    if surface_type not in builders:
        raise ValueError(f"不支持的表面类型 {surface_type!r}，可选: {', '.join(builders)}")
    if output_format not in ("vasp", "lammps-data"):
        raise ValueError(f"不支持的输出格式 {output_format!r}，可选: vasp, lammps-data")

    atoms = builders[surface_type](
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
        output_format="lammps-data",
    )

    # (110) 表面
    get_fcc_surface(
        symbol=symbol,
        a=a,
        surface_type="110",
        slab_size=(1, 1, 6),
        output_fn=file_110,
        output_format="lammps-data",
    )

    # (111) 表面
    get_fcc_surface(
        symbol=symbol,
        a=a,
        surface_type="111",
        slab_size=(1, 2, 6),
        output_fn=file_111,
        output_format="lammps-data",
    )

    surface_folder = Path("surfaces-ase")
    surface_folder.mkdir(exist_ok=True)
    for fn in (file_100, file_110, file_111):
        shutil.move(fn, surface_folder / fn)

    """
    from ase.data import atomic_numbers
    from ase.io import read

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
