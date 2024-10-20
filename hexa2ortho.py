"""六方转正交"""

from ase.build.supercells import make_supercell
from ase.io import read, write

# 六方转正交，reference: https://mp.weixin.qq.com/s/H66KxwJCLvr-hR-sosPB_w
# TODO: 该转变矩阵只适合夹角为 60 度的情况？


def hexa2ortho(input_fn: str, output_fn: str):

    T = [
        [1, 1, 0],
        [-1, 1, 0],
        [0, 0, 1],
    ]

    atoms_hexagonal = read(input_fn)

    atoms_orghogonal = make_supercell(
        prim=atoms_hexagonal,
        P=T,
    )

    write(
        output_fn,
        atoms_orghogonal,
        format="vasp",
        vasp5=True,
        sort=True,
        direct=True,
    )


if __name__ == "__main__":
    input_fn = "Zr_hexagonal.vasp"
    output_fn = "test.vasp"

    hexa2ortho(input_fn, output_fn)
