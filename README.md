# CMS Scripts

计算材料科学相关代码/脚本（脚本总数：100+）。

注意事项：

- 本仓库代码主要依赖 ASE、pymatgen、atomate、PyXtal、spglib 等 Python 包，以及 atomsk、vaspkit 等工具

- 使用脚本前，请先阅读源码及其中的注释！

- 部分脚本可解析命令行参数，使用 `python xxx.py -h` 或 `./xxx.py -h` 查看帮助信息

- 对于 ipynb 格式脚本文件，可使用 [nbviewer](https://nbviewer.org/) 在线查看

---

## 脚本内容

- `ase-usage/`: ASE 程序使用

```bash
ase_atoms.ipynb                      # atoms 模块使用
ase_spacegroup.ipynb                 # spacegroup 模块使用
ase_structure_generation.ipynb       # ase 构建晶体
ase_db.ipynb                         # db 模块使用
ase_data.ipynb                       # data 模块使用
ase_io.ipynb                         # io 模块使用
ase_visualize.ipynb                  # visualize 模块使用
ase_eos.ipynb                        # 拟合 EOS 方程并绘制曲线
ase_rdf.ipynb                        # 计算并绘制 RDF
ase_outcar.ipynb                     # 获取 OUTCAR 文件中数据（离子步构型、能量及受力）
ase_phase_diagram.ipynb              # 使用 ase phasediagram 模块绘制相图
surface_fcc.sh                       # FCC (100), (110), (111) 表面模型构建
surface_bcc.sh                       # BCC (100), (110), (111) 表面模型构建
read_mtp_cfg.py                      # 使用 ASE 读取 MTP 训练集中的 cfg 文件
```

---

- `pymatgen-usage/`: pymatgen 程序使用

```bash
pymatgen_periodic_table.ipynb        # periodic_table 模块使用；元素周期表 tui 版本绘制
pyamtgen_composition.ipynb           # pyamtgen composition 模块 使用

pymatgen_eos.ipynb                   # 拟合 EOS 方程并绘制曲线
pymatgen_unit.ipynb                  # unit 模块使用
pymatgen_symmetry.ipynb              # symmetry 模块使用；对称性分析
pymatgen_neighbor.ipynb              # 查看晶体结构中的原子近邻情况

pymatgen_slab.ipynb                  # surface 模块使用；构建表面 slab 模型
pymatgen_visualize_structure.ipynb   # 可视化晶体结构

pymatgen_vasp_help.ipynb             # pymatgen.io.vasp.help 模块使用
pymatgen_kpoints.ipynb               # pyamtgen.io.vasp.inputs 模块中的 Kpoints 类使用
pymatgen_Potcar.ipynb                # pyamtgen.io.vasp.inputs 模块中的 Potcar 类使用
pymatgen_vasprun.ipynb               # pymatgen.io.vasp.outputs 模块中的 Vasprun 类使用
pymatgen_outcar.ipynb                # pymatgen.io.vasp.outputs 模块中的 Outcar 类使用
pymatgen_oszicar.ipynb               # pymatgen.io.vasp.outputs 模块中的 Oszicar 类使用
pymatgen_vaspsets.ipynb              # pymatgen.io.vasp.sets 模块使用

pymatgen_dos.ipynb                   # pymatgen DOS 绘制
pymatgen_dos2.ipynb                  # pymatgen DOS 数据获取

pymatgen_elastic_properties.ipynb    # 使用 pymatgen ElasticTensor 类计算弹性性质

pymatgen_phase_diagram.ipynb         # phase diagram 模块使用；三元相图 convex hull 绘制

pymatgen_new_api.ipynb               # 新 api 模块使用（三元、四元相图 convex hull 绘制）
pymatgen_old_api.ipynb               # 旧 api 模块使用
```

---

- `VASP-scripts/`: VASP 相关脚本

```bash
VASP-Official-Tutorials-2017.pdf     # VASP 官方教程 2017 版；来自 https://github.com/tamaswells/VASP_script

check_force.py                       # 检查 OUTCAR 文件中离子步的原子受力收敛性；基于内置 re 模块
check_force_ase.py                   # 检查 OUTCAR 文件中离子步的原子受力收敛性；基于 ASE 程序
check_force_pymatgen.py              # 检查 OUTCAR 文件中离子步的原子受力收敛性；基于 pymatgen 程序
extract_force.sh                     # 提取离子步中的原子位置及受力，可指定原子、离子步步数
read_force.py                        # 解析每个目录下的 vasprun.xml 文件，提取受力并统计最大受力
extract_outcar.ipynb                 # 使用正则表达式提取 OUTCAR 文件中的数据（练习用）

CheckVaspDone                        # 检查 VASP 计算是否完成（可多个目录）
CheckOptConverged                    # 检查 VASP 弛豫计算是否收敛（可多个目录）
sigma.sh                             # 检查 entropy T*S 是否小于 1 meV/atom

potcar_pbe_compare.py                # 比较 VASP 和 pymatgen 推荐的常用元素 PBE 赝势 (VASP5.4.4)
get_potcar.py                        # 生成 VASP 和 pymatgen 推荐的赝势 POTCAR 文件
get_vasp_data.py                     # 获取 VASP 计算目录的输出数据（利用 atomate package）
get_vasp_data2.py                    # 获取 VASP 计算目录的输出数据（利用 pymatgen package）
get_vasp_data.sh                     # 获取多个 VASP 计算目录的输出数据
kgrid2kspacing.py                    # 从 KPOINTS 和 POSCAR 文件获取 K-spacing 数值
kspacing2kgrid.py                    # 从 K-spacing 数值和 POSCAR 获取 K-grid 数值
incar_generation.py                  # 根据计算类型生成 VASP INCAR 文件
get_kpts_atomate.py                  # 获取 atomate 中不同 workflows 的默认 kpts 设置
get_kpath.py                         # 根据构型生成 K-Path

plot_ev.py                           # EV（能量-体积）曲线绘制
eos_fit.py                           # Birch-Murnaghan EOS 拟合

hdf5.ipynb                           # 读取 vaspout.h5 格式数据文件内容（无法获取具体数值）

# BandStructure-DOS/ 目录
dos_bs_inputs.py                     # 生成 DOS 和 BandStructure 计算的输入文件
plot_bs_pymatgen.py                  # 使用 pymatgen 模块绘制能带结构
plot_bs_vaspkit.py                   # 使用 vaspkit 工具获取能带数据，绘制能带结构
get_dos_spd.py                       # 获取体系分态密度数据
plot_dos_spd.py                      # 使用 pymatgen 模块绘制体系分态密度
plot_dos_spd.py                      # 绘制体系分态密度（不使用 pymatgen 模块）
plot_dos_element_spd.py              # 使用 pymatgen 模块绘制元素分态密度
```

---

- `structure-scripts/`: 结构相关脚本

```bash
layers_count.py                      # 统计原子层数及每层原子数
interlayer_separations.py            # 统计原子层间距变化（适用于表面/界面模型弛豫前后的原子层间距变化）
identify_layer.py                    # 识别每个原子所在的原子层
fix_layer.py                         # 固定（特定）原子层 x/y/z 轴
wrap_pos.py                          # 将 VASP POSCAR 中分数坐标范围 wrap 在 0-1 之间
pos_diff.py                          # 比较构型弛豫前后原子坐标的变化
get_distance.py                      # 获取构型中原子对的最小和最大距离
posconv.py                           # 构型文件格式互相转换（基于 ASE，支持 ASE 大部分可识别的格式）
atat.py                              # 解析 ATAT 中的 str.out 文件（单个和枚举）的构型并转换为 ASE Atoms 对象
hexa2ortho.py                        # 将六方胞转换为正交胞
build_structure_spacegroup.ipynb     # PyXtal、ASE、pymatgen 通过空间群构建复杂结构
symmetry_info.py                     # 获取结构的对称性信息
sg_info.py                           # 获取结构的空间群信息
spglib_python.ipynb                  # 使用 spglib 的 Python API 获取对称性信息
spglib_julia.ipynb                   # 使用 spglib 的 Julia API 获取对称性信息
get_pearson_symbol.py                # 获取结构的 Pearson 符号
hcp_direction_index.py               # HCP 方向指数的三指数和四指数坐标转换
hcp_plane_index.py                   # HCP 面指数的三指数和四指数坐标转换
get_cn.py                            # 计算原子配位数


get_interface.py                     # 生成界面结构（待优化）

read_poscar.py                       # 读取 POSCAR 文件（练习用）
build_structure.py                   # 构建晶体结构
build_structure.ipynb                # 构建晶体结构（练习用）
```

---

- `NEP`: NEP 相关脚本

```bash
count_cfg.sh                         # 统计 cfg 构型文件帧数及总原子数
count_xyz.sh                         # 统计 extxyz 构型文件帧数及总原子数

json2extxyz.py                       # 将 json 构型及其数据文件转换为 extxyz 格式
son2db.py                            # 将 json 构型及其数据文件转换为 ASE db 格式
json2df.py                           # 将 json 构型及其数据文件转换为 Pandas DataFrame 格式

relax.py                             # 结构弛豫 ASE 实现
eos_cal_calorine.py                  # 使用 GPUMD & calorine 进行 EOS 计算
nep_pca_plot.py                      # NEP 未归一化 & 归一化 描述符 PCA 2 维绘制

plot_nep_loss.py                     # 绘制 NEP 训练 loss 演化与能量、力、virial/应力 RMSE 指标 DFT计算值 NEP 预测值对比图及其对应误差分布图
plt_nep_train.py                     # 绘制 NEP 训练 loss 演化与能量、力、virial/应力 RMSE 指标 DFT计算值 NEP 预测值对比图
```

---

- `ovito-usage/`: ovito Python 使用

```bash
get_point_defects.py                 # 获取轨迹文件中每帧构型的点缺陷（空位、间隙）数目
get_rdf_ovito.py                     # 计算轨迹文件平均 RDF
```

- `atomsk-usage/`: atomsk 程序使用

```bash
# 1-Surface/ 目录
surface_bcc.sh                       # BCC (100), (110), (111) 表面模型构建
surface_fcc.sh                       # FCC (100), (110), (111) 表面模型构建
surface_diamond.sh                   # Diamond (100), (110), (111) 表面模型构建

# 2-Stacking-Fault/ 目录
sf_fcc.sh                            # FCC 晶体的 ISF（本征堆垛层错）, ESF（非本征堆垛层错） 和 TWIN（孪晶）层错模型构建
twin_fcc.sh                          # FCC 孪晶模型构建

1-Al-polycrystal/                    # FCC Al 多晶模型构建
2-Fe-Cr-Ni-polycrystal/              # Fe-Cr-Ni 多晶模型构建
```

---

- `HPC/`: 超算用脚本

```bash
slurm_generation.py                  # 生成 VASP/LAMMPS/Python/Bash 任务 slurm 提交脚本
sruns                                # 根据超算类型申请计算节点
checkjob                             # 检查在队列中的 job 任务信息
print_help.sh                        # 定义 print_help 函数，供其他 Shell 脚本调用
```

---

- `atomate-usage/`: atomate 程序使用

```bash
static.py                            # 静态计算 workflow（默认参数，Si）
relax.py                             # 弛豫计算 workflow（默认参数，Si）
static_metal.py                      # 静态计算 workflow（修改 INCAR、KPOINTS 参数，以适用于金属体系）
relax_metal.py                       # 弛豫计算 workflow（修改 INCAR、KPOINTS 参数，以适用于金属体系）
atomate_basic.ipynb                  # atomate 基础使用
atomate_db.ipynb                     # 获取使用 atomate 计算并存储到数据库（Mongodb）中的数据
```

---

- `elastic/`: 弹性相关脚本

```bash
elastic_stability.py                 # Born 力学稳定性判据
elastic_property.ipynb               # 弹性性质计算
elastic_matrix.ipynb                 # 根据独立弹性常数和晶系生成弹性张量矩阵
```

---

- `sqsgen-usage/`: sqsgen 程序使用

---

- `vaspkit-usage/`: vaspkit 程序使用

```bash
Cu-bandstructure/                    # Cu 能带结构分析与绘制
Cu-DOS/                              # Cu DOS 分析与绘制
ELASTIC/                             # 从弹性张量文件中计算弹性性质（3D 和 2D）
```

---

- `atomkit-usage/`: atomkit 程序使用

```bash
Wyckoff.in                           # 查看对称性（Wyckoff）位置信息
```

---

- `image-process/`: 图片处理

```bash
pdf2img.py                           # pdf 转 png 格式图片
image_crop.py                        # 裁切图片多余空白
```

---

- `shell-scripts/`: Shell 脚本

```bash
pull.sh                              # 对个人常用的 Git 仓库进行批量 pull 操作
shell_set.ipynb                      # Shell set 命令测试
awk.ipynb                            # awk 命令使用
```

---

- `vasp6-mac-m1/`: 用 Apple M1 芯片编译 VASP6 过程中所涉及到的修改的源代码；查看修改细节，运行以下代码：

```bash
git diff cb92b24 cec4770
```

---

- `plots/`: 绘图示例

```bash
gamma_surface_plot.ipynb             # GSFE gamma surface 绘制
matplotlib_basic.ipynb               # 各种绘图示例
colors_template.py                   # 绘图配色
colors_template2.py                  # 绘图配色 2
corr_heatmap.ipynb                   # 相关性热图绘制
periodic_table_plot.ipynb            # 元素周期表绘制
misc_plot.ipynb                      # 其他绘图示例
```

---

- `misc/`: 其他脚本

```bash
clease.ipynb                         # clease 程序使用
dpdata.ipynb                         # dpdata 程序使用
pyxtal.ipynb                         # pyxtal 程序使用
pandas.ipynb                         # pandas 程序使用

periodic_table_info.py               # 打印元素周期表中元素的基本性质
coord_transform.ipynb                # 分数坐标与直角坐标互相转换
interatomic_potential_plot.ipynb     # 势函数示意图绘制
at2wt.py                             # 将 Ti-22Al-23Nb-1Mo-1Zr 格式化学式原子百分比转化成质量百分比
```
