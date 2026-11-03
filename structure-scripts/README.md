# structure-scripts

结构相关脚本

## 脚本列表

```bash
add_vacuum.py                        # 在指定轴（默认 z 轴）顶部、底部、顶部底部两端添加真空层（支持非正交胞，保持晶格角度不变）
layers_count.py                      # 统计原子层数及每层原子数
layers_identify.py                   # 识别原子层及其对应的原子
layers_fix.py                        # 固定（特定）原子层 x/y/z 轴
interlayer_separations.py            # 统计原子层间距变化（适用于表面/界面模型弛豫前后的原子层间距变化）

wrap_pos.py                          # 将 VASP POSCAR 中分数坐标范围 wrap 在 0-1 之间
pos_diff.py                          # 比较两个构型之间的原子坐标变化/差异（主要为弛豫/变形前后）
posconv.py                           # 构型文件格式互相转换（基于 ASE，支持 ASE 大部分可识别的格式）
supercell.py                         # 扩胞
atat.py                              # 解析 ATAT 中的 str.out 文件（单个和枚举）的构型并转换为 ASE Atoms 对象
hexa2ortho.py                        # 将六方胞转换为正交胞
build_structure_spacegroup.ipynb     # PyXtal、ASE、pymatgen 通过空间群构建复杂结构

symmetry_info.py                     # 获取结构的对称性信息
sg_info.py                           # 获取结构的空间群信息
get_pearson_symbol.py                # 获取结构的 Pearson 符号
spglib_python.ipynb                  # 使用 spglib 的 Python API 获取对称性信息
spglib_julia.ipynb                   # 使用 spglib 的 Julia API 获取对称性信息

hcp_direction_index.py               # HCP 方向指数的三指数和四指数坐标转换
hcp_plane_index.py                   # HCP 面指数的三指数和四指数坐标转换

get_cn.py                            # 计算原子配位数
get_nn.py                            # 获取 BCC/FCC/Diamond/HCP 晶体结构的最近邻距离
get_distance.py                      # 获取构型中原子对的最小和最大距离
get_orthogonal_c.py                  # 强制使 c 轴与 a、b 轴垂直（对于生成晶界、界面结构时有用）

element_data.py                      # 元素周期表常见元素数据

get_interface.py                     # 生成界面结构（待优化）

read_poscar.py                       # 读取 POSCAR 文件（练习用）
build_structure.py                   # 构建晶体结构
build_structure.ipynb                # 构建晶体结构（练习用）
```
