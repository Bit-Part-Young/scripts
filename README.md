# CMS Scripts

计算材料科学相关脚本。

注意事项：

- 本仓库代码主要依赖 ASE、pymatgen、atomate 等 Python 包

- 使用脚本前，请先阅读源码及其中的注释！

- 对于 ipynb 格式脚本文件，可使用 [nbviewer](https://nbviewer.org/) 在线查看

- 部分代码可解析命令行参数

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
surface_generation_fcc.sh            # FCC (100), (110), (111) 表面模型构建
read_mtp_cfg.py                      # 使用 ASE 读取 MTP 训练集中的 cfg 文件
```

---

- `pymatgen-usage/`: pymatgen 程序使用

```bash
pymatgen_periodic_table.ipynb        # periodic_table 模块使用；元素周期表 tui 版本绘制

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

check_force.py                       # 检查 OUTCAR 文件中的原子受力收敛性
check_force_ase.py                   # 检查 OUTCAR 文件中的原子受力收敛性；基于 ASE
read_force_pymatgen.ipynb            # 使用 pymatgen 读取 OUTCAR 文件中的原子位置与受力
extract_force.sh                     # 提取原子位置及受力，可指定原子、离子步步数
read_force.py                        # 解析每个目录下的 vasprun.xml 文件，提取受力并统计最大受力
extract_outcar.ipynb                 # 使用正则表达式提取 OUTCAR 文件中的数据（练习用）

CheckVaspDone                        # 检查 VASP 计算是否完成
CheckOptConverged                    # 检查 VASP 弛豫计算是否收敛（力 + 能量）
sigma.sh                             # 确定 entropy T*S 是否小于 1 meV/atom

potcar_pbe_compare.py                # 比较 VASP 和 pymatgen 推荐的常用元素 PBE 赝势 (VASP5.4.4)
get_potcar.py                        # 生成 VASP 和 pymatgen 推荐的赝势 POTCAR 文件
get_vasp_data.py                     # 获取 VASP 计算目录的数据

fit_ev.py                            # Birch-Murnaghan EOS 拟合

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
wrap_pos.py                          # 将 VASP POSCAR 中原子 Direct 坐标范围 wrap 在 0-1 之间
pos_diff.py                          # 比较构型弛豫前后原子坐标的变化
posconv.py                           # 构型文件格式互相转换（基于 ASE，支持 ASE 大部分可识别的格式）
atat.py                              # 解析 ATAT 中的 str.out 文件（单个和枚举）的构型并转换为 ASE Atoms 对象
hexa2ortho.py                        # 将六方胞转换为正交胞
build_structure_spacegroup.ipynb     # PyXtal、ASE、pymatgen 通过空间群构建复杂结构
crysinfo.py                          # 获取晶体结构对称性信息（WIP）
get_interface.py                     # 生成界面结构（待优化）
read_poscar.py                       # 读取 POSCAR 文件（练习用）
build_structure.ipynb                # 构建晶体结构（练习用）
```

---

- `atomsk-usage/`: atomsk 程序使用

```bash
1-Al-polycrystal/                    # FCC Al 多晶模型构建
2-Fe-Cr-Ni-polycrystal/              # Fe-Cr-Ni 多晶模型构建
surface_generation_fcc.sh            # FCC (100), (110), (111) 表面模型构建
surface_generation_bcc.sh            # BCC (100), (110), (111) 表面模型构建
twin_generation.sh                   # 孪晶模型构建
```

---

- `HPC/`: 超算用脚本

```bash
# 1-SiYuan、2-Pi 中的任务提交脚本
# 超算编译好的 LAMMPS 中安装的 package 很少，不推荐使用
vasp_local.slurm                     # 使用本地编译的 VASP
vasp.slurm                           # module load 加载超算编译好的 VASP
lammps_local.slurm                   # 使用本地编译的 LAMMPS
job_universal.slurm                  # 通用任务（Python、Shell 等）提交脚本

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

- 其他

```bash
clease.ipynb                         # clease 程序使用
dpdata.ipynb                         # dpdata 程序使用
pyxtal.ipynb                         # pyxtal 程序使用
pandas.ipynb                         # pandas 程序使用
```

---

- `shell-scripts/`: Shell 脚本

```bash
pull.sh                              # 对个人常用的 Git 仓库进行批量 pull 操作
shell_set.ipynb                      # Shell set 命令测试
awk.ipynb                            # awk 命令使用
```

---

- `misc/`: 其他脚本

```bash
coord_transform.ipynb                # 分数坐标与直角坐标互相转换
index_transform.ipynb                # 三、四指数坐标转换
interatomic_potential_schematic.ipynb  # 势函数示意图绘制
at2wt.ipynb                          # 将 Ti-22Al-23Nb-1Mo-1Zr 格式化学式原子百分比转化成质量百分比
```

---

- `vasp6-mac-m1/`: 用 Apple M1 芯片编译 VASP6 过程中所涉及到的修改的源代码；查看修改细节，运行以下代码：

```bash
git diff cb92b24 cec4770
```

---

- `plots/`: 绘图示例

```bash
plot_ev.py                           # EV（能量-体积）曲线绘制
matplotlib_basic.ipynb               # 各种绘图示例
colors.ipynb                         # 绘图配色
corr_heatmap.ipynb                   # 相关性热图绘制
eriodic_table_plot.ipynb             # 元素周期表绘制
misc_plot.ipynb                      # 其他绘图示例
```

---

- `deprecated/`: 弃用脚本

```bash
xsd2vasp.py                          # xsd POSCAR 构型文件格式互相转换
at2wt.py                             # 将 Ti-22Al-23Nb-1Mo-1Zr 格式化学式原子百分比转化成质量百分比
```
