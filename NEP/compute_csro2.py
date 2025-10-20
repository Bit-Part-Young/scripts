#!/usr/bin/env python3

"""
计算 CSRO / WCP 参数
reference: https://mp.weixin.qq.com/s/Yc4zXrPpRZp5yLX7Weksgg

- [ ] 代码待优化
"""

import warnings
from math import sqrt

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from ovito.data import CutoffNeighborFinder as CNF
from ovito.data import DataCollection
from ovito.io import import_file

# lattice_type = "fcc"  # 定义模型的晶格类型
# lattice = 3.6  # 定义晶格长度, 单位和dump文件中的长度单位一致
fsize = 24  # 绘制皮尔逊图的时候, 每个格子里面的数字字体和其他地方的字符串大小
pad = 0.03  # 刻度字符串距离边线的距离
lattice_type = "bcc"  # 定义模型的晶格类型
lattice = 3.289  # 定义晶格长度, 单位和dump文件中的长度单位一致


def find_cutoff(lat_type, a0):
    """
    此函数根据传入的晶格类型和晶格常数计算对应的截断半径
    输入参数:
    lat_type: 字符串类型, 目前只支持fcc和bcc的晶格类型
    a0:       浮点数, 对应当前模型的晶格常数
    输出参数:
    cutoff:   浮点数, 根据传入的晶格类型计算截断半径
    """
    # 参考资料:
    # https://docs.lammps.org/compute_cna_atom.html
    # https://linkinghub.elsevier.com/retrieve/pii/0927025694901090
    if lat_type == "fcc":
        cutoff = 0.25 * (sqrt(2) + 2.0) * a0
        # BCC 的 cutoff 有问题
    elif lat_type == "bcc":
        cutoff = 0.5 * (sqrt(2) + 1.0) * a0
    else:
        raise ValueError(
            "Unrecognized lattice type, please input correct lattice type."
        )

    return cutoff


def fprintf(WCP):
    """
    将计算出来的WCP格式化输出, 即保留2位小数, 输入参数:
    WCP: 计算完成之后的WCP系数矩阵
    write_to_file: 逻辑变量, 如果设置为真, 那么就会将结果输出到文件里面, 文件名称为result.txt
    """
    np.set_printoptions(
        formatter={"float": "{: 0.2f}".format}
    )  # 设置保留2位小数, 并且对很大/小的数不使用科学计数法
    num = len(WCP)  # 获取模型有几种元素
    transform_WCP = np.zeros((num, num), dtype=np.double)
    for i in range(num - 1, -1, -1):  # 循环输出WCP系数矩阵, 并配合下标进行输出
        temp = WCP[i, :]
        print(i + 1, temp)
        transform_WCP[abs(i - (num - 1)), :] = temp
    for i in range(1, num + 1):  # 输出下标, 这样更方便阅读
        print("   {}".format(i), end="   ")
    return transform_WCP


def plot_WCP(WCP):
    """
    本函数用于将传入的WCP系数矩阵绘制成对应的皮尔逊相关性热图
    """
    num = len(WCP)
    # 计算x轴和y轴刻度的字符串, 这里注意由于dump文件里面不一定会指出元素的名称
    # 所以这里统一用1号原子, 2号原子进行标识
    xlabel_ticks = list(range(0, len(WCP) + 1))
    ylabel_ticks = list(range(len(WCP) + 1, 0, -1))
    # 设置画布
    fig, ax = plt.subplots()
    # 绘制热图, 这里的cmap就是指定采用哪种预先内置的colorbar, 同学们有需要可以到以下网址查询
    # https://matplotlib.org/1.5.3/gallery.html , 同学们注意, 这里我提前设置了color bar的上下限为1.0
    # 和-2.0, 一般来说大多数模型应该都是位于这个区间里面, 同学们如果有特殊需求需要手动修改代码
    im = ax.imshow(
        WCP, cmap=plt.cm.viridis, interpolation="nearest", vmin=-2.0, vmax=1.0
    )
    # 这里我们设置右边的color bar的刻度数字都保留至少一位小数
    formatter = mpl.ticker.StrMethodFormatter("{x:.1f}")
    # 设置x和y轴的具体刻度字符串, 并设置字符串距离边框的距离
    ax.yaxis.set_major_locator(MaxNLocator(num + 1))
    ax.xaxis.set_major_locator(MaxNLocator(num + 1))
    ax.set_xticklabels(xlabel_ticks, fontsize=fsize, position=(0, -pad))
    ax.set_yticklabels(ylabel_ticks, fontsize=fsize, position=(-pad, 0))
    # 绘制colorbar , 并调用之前我们设置的格式化变量formatter
    temp = fig.colorbar(im, format=formatter, ax=ax)
    # 调整colorbar刻度线长度为5, 且字体大小为变量fsize的大小, 如果不想要刻度线, 可以将length设为0即可
    temp.ax.tick_params(length=5, labelsize=fsize)
    # 这里我们取消colorbar的外框线的可视性, 同学们根据自己需要进行自定义
    temp.outline.set_visible(False)
    # 这里我强制定义右边的colorbar只有3个刻度, 即上下限-2.0和1.0, 以及中间的0.0,
    # 同学们可以根据自己喜欢进行定制
    temp.set_ticks([-2.0, 0.0, 1.0])
    # 强制更新刚刚在colorbar上面所做的修改
    temp.update_ticks()
    # 在colorbar的旁边加一行字, 但是感觉效果不太好, 就注释掉了, 同学们感兴趣可以自行调整
    # temp.set_label("Warren-Cowley parameters", fontsize = fsize)

    # 设置标题
    ax.set_title("Warren-Cowley parameters(WCP) for current model")
    # 关掉刻度线的显示
    ax.tick_params(bottom=False, top=False, left=False, right=False)
    # 关掉边框的显示
    for i in ["top", "right", "bottom", "left"]:
        ax.spines[i].set_visible(False)
    # 双层循环, 在每个格子里面写入具体的WCP系数的值, 颜色通过color关键字定义, 这里默认是白色
    # 然后居中对齐, 字体大小通过变量fsize进行定义
    for i in range(len(WCP)):
        for j in range(len(WCP)):
            temp_str = "{:.2f}".format(
                WCP[i, j]
            )  # 将WCP矩阵每个系数保留2位小数的字符串
            text = ax.text(
                j, i, temp_str, fontsize=fsize, ha="center", va="center", color="w"
            )

    # fig.tight_layout()
    fig.savefig("WCP.png")

    # 显示图片
    plt.show()


def wcp_modifier(input: DataCollection):
    cutoff = find_cutoff(lattice_type, lattice)  # 计算截断半径
    finder = CNF(cutoff, input)  # 以计算出来的截断半径构建邻居列表
    type_array = input.particles["Particle Type"].array  # 获取类型数组
    ntypes = max(type_array)  # 计算模型中共有几种原子
    ratio = np.zeros(ntypes, dtype=np.double)  # 存储每种元素的百分比
    counts = np.zeros(
        ntypes, dtype=np.double
    )  # 存储每种元素在指定截断半径以内一共有多少个原子
    nparticles = input.particles.count  # 模型中总共有多少个原子
    WCP = np.zeros((ntypes, ntypes), dtype=np.double)  # 存储模型最终的WCP系数的矩阵

    # 调用循环计算每种元素的具体百分比
    for i in range(ntypes):
        ratio[i] = np.sum(type_array == (i + 1)) / nparticles

    # 调用循环计算每个元素的近邻原子数, 第一层循环遍历每个原子, 第二层循环遍历第i个原子的所有邻居原子
    for i in range(nparticles):
        cen_type = type_array[i]  # 获取中心原子的原子类型
        for neigh in finder.find(i):
            neigh_type = type_array[neigh.index]  # 获取中心原子的邻居原子的原子类型
            counts[cen_type - 1] += 1  # 中心原子的邻居原子数目加一
            WCP[
                cen_type - 1, neigh_type - 1
            ] += 1.0  # 统计中心原子的邻居原子中neigh_type的数量

    # 根据公式进行最后计算并根据对称性进行转置
    for i in range(ntypes):
        for j in range(i, ntypes):
            WCP[i, j] = 1.0 - WCP[i, j] / (ratio[j] * counts[i])
            if i != j:
                WCP[j, i] = WCP[i, j]
    WCP = fprintf(WCP)  # 格式化输出WCP矩阵
    # plot_WCP(WCP)  # 根据WCP矩阵绘制相关性热图


if __name__ == "__main__":

    # pipeline = import_file("CoCuFeNiPd-2M.txt")
    pipeline = import_file("dump.xyz")
    for frame in range(pipeline.num_frames):
        if frame > 10:
            continue
        data = pipeline.compute(frame)
        wcp_modifier(data)
        print(f"Frame {frame} processed.")
