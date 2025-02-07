#!/usr/bin/env python3

"""计算原子配位数"""

from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.core.structure import Structure


def get_coordination_number(structure_fn: str = "POSCAR", cutoff: float = 3.0):
    """计算原子配位数"""

    structure = Structure.from_file(structure_fn)
    frac_coords = structure.frac_coords.round(4).tolist()

    for i in range(len(structure)):
        voronoi = VoronoiNN(cutoff=cutoff)
        coordination_number = voronoi.get_cn(structure, n=i)

        info = [i, structure[i].specie.name, frac_coords[i], coordination_number]

        print("Atom {0} {1} {2} CN: {3}".format(*info))


if __name__ == "__main__":
    get_coordination_number()
