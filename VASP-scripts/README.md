# VASP Scripts

- `OUTCAR_test`: https://github.com/K4ys4r/VASP_Scripts/blob/master/Extract_Forces/OUTCAR_test

## 脚本列表

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

plot_energy_force.sh                 # 获取 VASP 弛豫过程中的能量、原子受力（考虑原子位置方向固定情况，计算受力的模长）信息并绘制演化图
plot_energy_force2.sh                # 获取 VASP 弛豫过程中的能量、原子受力（不考虑原子位置方向固定情况，不计算受力的模长）信息并绘制演化图
plot_rdf_AIMD.sh                     # 获取 AIMD 的 RDF 数据并绘制
plot_temperature_energy_AIMD.sh      # 获取 AIMD 的温度和能量数据并绘制
plot_ev.py                           # EV（能量-体积）曲线绘制

get_psp.py                           # 生成 VASP 和 pymatgen 推荐的 PBE 赝势 POTCAR 文件
get_psp2.py                          # 生成 VASP 和 pymatgen 推荐的 PBE 赝势 POTCAR 文件（不依赖 pymatgen，速度更快）
psp_comparison.py                    # 比较 VASP 和 pymatgen 推荐的常用元素 PBE 赝势 (VASP5.4.4)
get_vasp_data.py                     # 获取 VASP 计算目录的输出数据（利用 atomate package）
get_vasp_data2.py                    # 直接解析 VASP 能量、离子步数及实际耗时（无需 pymatgen）；支持多个目录、. 和 '*'
get_opt_data.sh                      # 获取当前 目录/子目录 下的 VASP 弛豫计算数据
get_scf_data.sh                      # 获取当前 目录/子目录 下的 VASP 静态计算数据
get_vasp_running_job.sh              # 获取正在运行的 VASP 作业信息
vasp_timecost.sh                     # 统计 VASP 计算目录的离子步步数、计算耗时、能量信息
kgrid2kspacing.py                    # 从 KPOINTS 和 POSCAR 文件获取 K-spacing 数值
kspacing2kgrid.py                    # 从 K-spacing 数值和 POSCAR 获取 K-grid 数值
incar_generation.py                  # 根据计算类型生成 VASP INCAR 文件
get_kpts_atomate.py                  # 获取 atomate 中不同 workflows 的默认 kpts 设置
get_kpath.py                         # 根据构型生成 K-Path
eos_fit.py                           # Birch-Murnaghan EOS 拟合
get_surface_energy.py                # 计算表面能
pde_cal.py                           # 点缺陷形成能计算

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
