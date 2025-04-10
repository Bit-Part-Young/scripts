"""
自动寻找间隙原子位置

BCC 找到了 四面体间隙
FCC 找到了 四、八面体间隙
HCP 找到了 四、八面体间隙（四面体间隙 z 轴分数坐标不太一致，应为 3/4、1/4，对于原胞而言）

reference: https://github.com/bracerino/Automatically-find-interstitial-sites/blob/main/calc_int_sites.py
"""

import numpy as np
from pymatgen.analysis.defects.core import Interstitial
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp import Poscar

np.set_printoptions(precision=5, suppress=True)


interstitial_element = "N"
number_interstitials_insert = 1
# 0 表示考虑所有类型的间隙原子，1 表示只考虑第一种间隙原子，2 表示只考虑第二种间隙原子
interstitial_type_selected = 0


structure = Poscar.from_file("POSCAR").structure
generator = VoronoiInterstitialGenerator(
    # clustering_tol=0.5,
    # clustering_tol=1.0,
    clustering_tol=0.75,
    # min_dist=0.9,
    min_dist=0.5,
)


# 获取所有类型的间隙原子及对应的等同位置
interstitial_site_total_list = []
unique_interstitial_list = []
unique_multiplicity_list = []
interstitial_type = 0
interstitial_site_dict = {}
# Element 'H' is here only to find the available sites in order to prevent error with oxidation states for some elements like noble gases
for interstitial in generator.generate(structure, "H"):
    interstitial: Interstitial
    interstitial_site_dict[interstitial_type] = []

    interstitial_site = interstitial.site.frac_coords
    multiplicity = interstitial.multiplicity
    print(f"Unique interstitial site: {interstitial_site}, multiplicity: {multiplicity}")
    print(f"\nThe equivalent positions:")

    for index, equivalent_site in enumerate(interstitial.equivalent_sites, start=1):
        print(f"Position {index:02d}: {equivalent_site.frac_coords}")
        interstitial_site_total_list.append(equivalent_site.frac_coords)
        interstitial_site_dict[interstitial_type].append(equivalent_site.frac_coords)

    print(f"\n---------------------------------------------------------\n")

    unique_interstitial_list.append(interstitial_site)
    unique_multiplicity_list.append(multiplicity)

    interstitial_type = interstitial_type + 1

print(f"There are total of {len(unique_interstitial_list)} unique interstitial sites:")
for equivalent_site, multiplicity in zip(unique_interstitial_list, unique_multiplicity_list):
    print(f"site: {equivalent_site}, multiplicity: {multiplicity}")


if interstitial_type_selected == 0:
    frac_coords_use = interstitial_site_total_list
else:
    frac_coords_use = interstitial_site_dict[interstitial_type_selected - 1]


def wrap_coordinates(frac_coords):
    """Wrap fractional coordinates into the range [0, 1)."""

    frac_coords = np.array(frac_coords)  # Ensure input is a NumPy array
    return frac_coords % 1


def distance_matrix_cal_pbc(frac_coords: np.ndarray) -> np.ndarray:
    """计算原子间距（考虑 PBC）"""

    n = len(frac_coords)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            delta = frac_coords[i] - frac_coords[j]
            delta -= np.round(delta)
            dist_matrix[i, j] = dist_matrix[j, i] = np.linalg.norm(delta)

    return dist_matrix


def select_spaced_points(
    frac_coords,
    n_points=5,
    mode="farthest",
    target_value=0.5,
):
    """选择 n 个等间距的点"""

    frac_coords = wrap_coordinates(frac_coords)
    dist_matrix = distance_matrix_cal_pbc(frac_coords)

    selected_indices = [0]
    for _ in range(1, n_points):
        remaining_indices = [i for i in range(len(frac_coords)) if i not in selected_indices]

        if mode == "farthest":
            next_index = max(
                remaining_indices,
                key=lambda i: min(dist_matrix[i, j] for j in selected_indices),
            )
        elif mode == "nearest":
            next_index = min(
                remaining_indices,
                key=lambda i: min(dist_matrix[i, j] for j in selected_indices),
            )
        elif mode == "moderate":
            next_index = min(
                remaining_indices,
                key=lambda i: abs(
                    sum(dist_matrix[i, j] for j in selected_indices) / len(selected_indices)
                    - target_value
                ),
            )
        else:
            raise ValueError("Invalid mode. Choose from 'nearest', 'farthest', or 'moderate'.")

        selected_indices.append(next_index)

    return frac_coords[selected_indices].tolist()


def insert_interstitials(frac_coords, n_points=5, mode="farthest"):
    selected_points = select_spaced_points(
        frac_coords,
        n_points=n_points,
        mode=mode,
    )

    structure_copy = structure.copy()
    for point in selected_points:
        structure_copy.append(
            species=Element(interstitial_element),
            coords=point,
            coords_are_cartesian=False,
        )

    poscar = Poscar(structure_copy)
    poscar.write_file(f"interstitials_{mode}.vasp")
