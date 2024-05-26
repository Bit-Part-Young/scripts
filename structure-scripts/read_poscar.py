import numpy as np


def read_POSCAR(filename):
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()

        comment = lines[0].strip()
        scale = float(lines[1].strip())
        lattice = np.array(
            [list(map(float, x.strip().split())) for x in lines[2:5]]
        )
        symbols = lines[5].strip().split()
        natoms_list = list(map(int, lines[6].strip().split()))

        if lines[7].strip() in ["s", "S", "Selective dynamics"]:
            coords = lines[9:]
            coord_type = lines[8].strip()
        else:
            coord_type = lines[7].strip()
            coords = lines[8:]

        coord = np.array(
            [list(map(float, x.strip().split()[:3])) for x in coords]
        )
        selective_dynamics = [x.strip().split()[3:] for x in coords]
        selective_dynamics = np.array(
            [[1 if i == "T" else 0 for i in x] for x in selective_dynamics]
            # 将 selective_dynamics 转换为 bool 类型
            # [[True if i == "T" else False for i in x] for x in selective_dynamics]
        )

    return (
        comment,
        scale,
        lattice,
        symbols,
        natoms_list,
        coord_type,
        coord,
        selective_dynamics,
    )


if __name__ == "__main__":
    fn = "POSCAR"
    poscar = read_POSCAR(filename=fn)
    print(poscar)
