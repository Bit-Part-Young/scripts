"""六方胞转正交胞"""

from ase.build.supercells import make_supercell
from ase.io import read, write

# 六方转正交，reference: https://mp.weixin.qq.com/s/H66KxwJCLvr-hR-sosPB_w

# TODO: 如何将转换后的 正交胞 lattice vector 转换成只有 xx yy zz 分量的格式


def hexa2ortho(input_fn: str, output_fn: str):
    """六方胞转正交胞"""

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
    input_fn = "hexagonal.vasp"
    output_fn = "orthogonal.vasp"

    hexa2ortho(
        input_fn=input_fn,
        output_fn=output_fn,
    )
