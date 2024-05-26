import os

import numpy as np
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder

# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as Analyzer
from pymatgen.analysis.interfaces.zsl import ZSLGenerator
from pymatgen.core.structure import Structure
from pymatgen.core.surface import get_symmetrically_distinct_miller_indices


def Get_Interface(
    prefix,
    subs,  # 基体相Nb
    film,  # 析出相Nb5Si3
    subs_miller,
    film_miller,
    saved_poscar_path,
    termination=None,
    max_area_ratio_tol=0.06,  # 四个max相关参数 生成匹配heterostructural interfaces（共格界面）需用到的最大容差参数
    max_area=100.0,
    max_length_tol=0.03,
    max_angle_tol=0.01,
    substrate_thickness=3,  # 构建基体相Nb slab 含原子层/layer 的最小尺寸  SlabGenerator类中的min_slab_size参数 默认值为1埃
    film_thickness=3,  # 构建析出相Nb5Si3 slab 含原子层/layer 的最小尺寸  SlabGenerator类中的min_slab_size参数 默认值为1埃
    gap=2,  # 两个相/slab 之间的距离 Interface.from_slabs()中的gap参数
    vacuum_over_film=15,  # 析出相上方的真空层厚度 CoherentInterfaceBuilder类中get_interfaces()方法中的vacuum_over_film参数
    frozen_thickness=5,
):

    zslgen = ZSLGenerator(
        max_area_ratio_tol=max_area_ratio_tol,
        max_area=max_area,
        max_length_tol=max_length_tol,
        max_angle_tol=max_angle_tol,
        bidirectional=False,
    )

    builder = CoherentInterfaceBuilder(
        substrate_structure=subs,
        film_structure=film,
        substrate_miller=subs_miller,
        film_miller=film_miller,
        zslgen=zslgen,
    )

    count = 0

    # if termination == None:
    if not termination:
        termination = builder.terminations[0]  # 析出相的表面终端
    # termination = builder.terminations[0]  # 析出相的表面终端

    for interface in builder.get_interfaces(
        termination=termination,
        substrate_thickness=substrate_thickness,
        film_thickness=film_thickness,
        gap=gap,
        vacuum_over_film=vacuum_over_film,
        in_layers=True,
    ):  # 使用层数作为厚度单位

        count = count + 1

        # 将未进行原子固定处理的界面结构的基体相和析出相slab分别保存成poscar文件
        fname_substrate = f"{saved_poscar_path}/{prefix}_{count}_subs.poscar"
        interface.substrate.to(fmt="poscar", filename=fname_substrate)

        fname_film = f"{saved_poscar_path}/{prefix}_{count}_film.poscar"
        interface.film.to(fmt="poscar", filename=fname_film)

        # 将未进行原子固定处理的界面结构保存成poscar文件
        fname_interface_unfixed = f"{saved_poscar_path}/{prefix}_{count}_unfixed.poscar"
        interface.to(fmt="poscar", filename=fname_interface_unfixed)

        # num_sites 界面的原子数  给界面所有原子设置[1 1 1]的坐标形式，方便后面的numpy计算 以及表示所有原子都是非固定的 True
        sd = np.ones((interface.num_sites, 3))
        # cart_coords 界面所有原子的笛卡尔坐标 axis=0 [2]表示对坐标在列方向上且是z轴上的坐标进行操作
        zmin = np.min(interface.cart_coords, axis=0)[2] + frozen_thickness
        for i, site in enumerate(interface.sites):  # 界面原子位点的iterator
            if site.coords[2] < zmin:
                # 固定的原子[1 1 1]坐标变成[0 0 0] 相当于对远离界面区域的原子进行固定 基体相slab
                sd[i][0] = sd[i][1] = sd[i][2] = 0
        # structure.add_site_property("selective_dynamics", [[True, True, True], [False, False, False]])
        interface.add_site_property("selective_dynamics", sd)

        natoms = interface.num_sites
        nfixed = natoms - np.sum(sd) / 3  # np.sum(sd) / 3 为非固定的原子数
        ntmsub = interface.substrate.num_sites  # 基体相Nb slab的原子数
        ntmflm = interface.film.num_sites  # 析出相Nb5Si3 slab的原子数
        subthick = "{:.3f}".format(
            np.max(interface.substrate.cart_coords, axis=0)[2]
            - np.min(interface.substrate.cart_coords, axis=0)[2]
        )
        flmthick = "{:.3f}".format(
            np.max(interface.film.cart_coords, axis=0)[2]
            - np.min(interface.film.cart_coords, axis=0)[2]
        )

        zmax = np.max(interface.cart_coords, axis=0)[2] - frozen_thickness
        for i, site in enumerate(interface.sites):
            if site.coords[2] > zmax:
                # 相当于对远离界面区域的原子进行固定 析出相slab
                sd[i][0] = sd[i][1] = sd[i][2] = 0
        interface.add_site_property("selective_dynamics", sd)

        # 将进行原子固定处理后的界面结构保存成poscar文件
        fname_interface = f"{saved_poscar_path}/{prefix}_{count}_fixed.poscar"
        interface.to(fmt="poscar", filename=fname_interface)

        # 'substrate_sl_vectors' 和 'film_sl_vectors' 为基体相和析出相slab超点阵/超胞的x y方向矢量
        # 这两个变量在ZSLMatch类中涉及 ZSLGenerator类有调用ZSLMatch
        su_sl = subs.lattice.get_fractional_coords(
            interface.interface_properties["substrate_sl_vectors"][0]
        )
        sv_sl = subs.lattice.get_fractional_coords(
            interface.interface_properties["substrate_sl_vectors"][1]
        )
        fu_sl = film.lattice.get_fractional_coords(
            interface.interface_properties["film_sl_vectors"][0]
        )
        fv_sl = film.lattice.get_fractional_coords(
            interface.interface_properties["film_sl_vectors"][1]
        )
        su_sl = "[{0:.1f},{1:.1f},{2:.1f}]".format(*su_sl)
        sv_sl = "[{0:.1f},{1:.1f},{2:.1f}]".format(*sv_sl)
        fu_sl = "[{0:.1f},{1:.1f},{2:.1f}]".format(*fu_sl)
        fv_sl = "[{0:.1f},{1:.1f},{2:.1f}]".format(*fv_sl)

        area = "{0:.4f}".format(interface.lattice.volume / interface.lattice.c)
        # 'strain'和'von_mises_strain'在  CoherentInterfaceBuilder类中的get_interfaces()中涉及
        voigt = interface.interface_properties["strain"].voigt
        strain = np.sqrt(np.sum(voigt[0:4] * voigt[0:4]) / 3)
        vmises = interface.interface_properties["von_mises_strain"]
        strain = "{0:.4f}".format(strain)
        vmises = "{0:.4f}".format(vmises)

        lattice_substrate = np.array(list(interface.substrate.lattice.abc))
        # lattice_substrate = np.array([interface.substrate.lattice[i][i] for i in range(3)])
        lattice_film = np.array(list(interface.film.lattice.abc))
        lattice_interface = np.array(list(interface.lattice.abc))

        misfit_lattice_sub_in = (
            (lattice_interface - lattice_substrate) / lattice_substrate
        )[0:2]
        misfit_lattice_film_in = ((lattice_interface - lattice_film) / lattice_film)[
            0:2
        ]

        # print(termination)
        print(
            f"{count}-{prefix} {lattice_interface} {lattice_substrate} {lattice_film} {misfit_lattice_sub_in} {misfit_lattice_film_in}"
        )
        print(
            f"{natoms}({ntmsub}+{ntmflm})-nfixed: {nfixed:.0f} {subthick} {flmthick} {su_sl} {sv_sl} {fu_sl} {fv_sl} {area} {vmises} {strain}"
        )
        print("\n")

    print(f"{count} interfaces is generated.")


def get_miller_index(subs, film, max_index):
    sms = get_symmetrically_distinct_miller_indices(subs, max_index)
    fms = get_symmetrically_distinct_miller_indices(film, max_index)

    return sms, fms


# 基体相与析出相的slab表面index已知 可不使用该函数
def Find_Interfaces(
    prefix,
    subs,
    film,
    subs_miller,
    film_miller,
    saved_poscar_path,
    max_area_ratio_tol=0.06,
    max_area=100.0,
    max_length_tol=0.03,
    max_angle_tol=0.01,
    substrate_thickness=4.0,
    film_thickness=1,
    gap=2,
    vacuum_over_film=0,
    frozen_thickness=5,
):
    zslgen = ZSLGenerator(
        max_area_ratio_tol=max_area_ratio_tol,
        max_area=max_area,
        max_length_tol=max_length_tol,
        max_angle_tol=max_angle_tol,
        bidirectional=False,
    )

    builder = CoherentInterfaceBuilder(
        substrate_structure=subs,
        film_structure=film,
        substrate_miller=subs_miller,
        film_miller=film_miller,
        zslgen=zslgen,
    )
    i = 0
    for term in builder.terminations:
        print(term)
        i += 1
        Get_Interface(
            prefix
            + f"_{''.join(map(str, subs_miller))}_{''.join(map(str, film_miller))}_termination_{i}",
            subs,
            film,
            subs_miller,
            film_miller,
            saved_poscar_path,
            termination=term,
            max_area_ratio_tol=max_area_ratio_tol,
            max_area=max_area,
            max_length_tol=max_length_tol,
            max_angle_tol=max_angle_tol,
            substrate_thickness=substrate_thickness,
            film_thickness=film_thickness,
            gap=gap,
            vacuum_over_film=vacuum_over_film,
            frozen_thickness=frozen_thickness,
        )


if __name__ == "__main__":
    prefix = "Nb-Nb5Si3"
    subs = Structure.from_file("POSCAR_Nb")
    film = Structure.from_file("POSCAR_Nb5Si3-alpha")
    subs_miller = (0, 0, 1)
    film_miller = (0, 0, 1)
    # substrate_thickness = 2
    # film_thickness = 2

    # sms, fms = get_miller_index(subs, film, max_index=1)
    # print(sms, fms)
    # subs_miller_index [(1, 1, 1), (1, 1, 0), (1, 0, 0)]
    # film_miller_index [(1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0), (0, 0, 1)]

    saved_poscar_path = "interface-poscar-6"
    os.path.exists(saved_poscar_path) or os.makedirs(saved_poscar_path)
    # Get_Interface(prefix, subs, film, subs_miller, film_miller,
    #               saved_poscar_path, substrate_thickness, film_thickness)
    # 界面结构的超胞大小是2x2x2Nb超胞的大小

    Find_Interfaces(prefix, subs, film, subs_miller, film_miller, saved_poscar_path)

    print("\nwork is done.\n")
