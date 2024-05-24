# CMS Scripts

计算材料科学相关脚本。

---

## 脚本内容

>对于 ipynb 格式文件，可使用 [nbviewer](https://nbviewer.org/) 在线查看。

- `ase-usage/`: ASE 程序使用
  - `ase_atoms.ipynb`: atoms 模块使用
  - `ase_db.ipynb`: db 模块使用
  - `ase_io.ipynb`: io 模块使用
  - `ase_visualize.ipynb`: visualize 模块使用
  - `ase_crystal.ipynb`: 构建晶体
  - `ase_eos.ipynb`: 拟合 eos 方程并绘制曲线
  - `ase_rdf.ipynb`: 计算并绘制 rdf
  - `ase_outcar.ipynb`: 获取 OUTCAR 文件中数据（离子步构型、能量及受力）

---

- `pymatgen-usage/`: pymatgen 程序使用
  - `pymatgen_eos.ipynb`: 拟合 eos 方程并绘制曲线
  - `pymatgen_periodic_table.ipynb`: periodic_table 模块使用；元素周期表 tui 版本绘制
  - `pymatgen_symmetry.ipynb`: symmetry 模块使用；对称性分析
  - `pymatgen_phase_diagram.ipynb`: phase diagram 模块使用；三元相图 convex hull 绘制
  - `pymatgen_api.ipynb`: 新旧 api 模块使用
  - `pymatgen_neighbor.ipynb`: 查看晶体结构中的原子近邻情况
  - `pymatgen_slab.ipynb`: surface 模块使用；构建表面 slab 模型
  - `pymatgen_vasp_help.ipynb`: `pymatgen.io.vasp.help` 模块使用
  - `pymatgen_vasp_outcar.ipynb`: `pymatgen.io.vasp.outputs` 模块中的 `Outcar` 类使用
  - `pymatgen_vasp_oszicar.ipynb`: `pymatgen.io.vasp.outputs` 模块中的 `Oszicar` 类使用
  - `pymatgen_elastic_properties.ipynb`: 使用 pymatgen ElasticTensor 类计算弹性性质

---

- `atomate-usage/`: atomate 程序使用
  - `atomate_db.ipynb`: 获取使用 atomate 计算并存储到数据库（Mongodb）中的数据

---

- `VASP-scripts`: VASP 相关脚本
  - `check_force.py`: 检查 OUTCAR 文件中的原子受力收敛性
  - `check_force_ase.py`: 检查 OUTCAR 文件中的原子受力收敛性；基于 ASE
  - `uniform_direct_coords.py`: 将 VASP POSCAR 中原子坐标分数范围限制在 0-1 之间
  - `read_force_pymatgen.ipynb`: 使用 pymatgen 读取 OUTCAR 文件中的原子位置与受力
  - `extract_force.sh`: 提取原子位置及受力，可指定原子、离子步步数
  - `read_force.py`: 解析每个目录下的 vasprun.xml 文件，提取受力并统计最大受力

---

- `shell-scripts/`: shell 脚本
  - `shell_set.ipynb`: shell set 命令测试

---

- `image-process/`: 图片处理
  - `pdf2img.py`: pdf 转 png 格式图片

---

- 其他
  - `posconv.py`: 构型文件格式互相转换；基于 ASE
  - `atat.py`: 解析 ATAT 中的 str.out 格式构型文件；基于 ASE
  - `clease.ipynb`: clease 程序使用
  - `dpdata.ipynb`: dpdata 程序使用
  - `pyxtal.ipynb`: pyxtal 程序使用

- `misc/`: 其他脚本
  - `at2wt.py`: 将 Ti-22Al-23Nb-1Mo-1Zr 格式化学式原子百分比转化成质量百分比
  - `coord_transform.ipynb`: 分数坐标与直角坐标互相转换
  - `index_transform.ipynb`: 三、四指数坐标转换
  - `interatomic_potential_schematic.ipynb`: 势函数示意图绘制

---

- `vasp6-mac-m1/`: 用 Apple M1 芯片编译 VASP6 过程中所涉及到的修改的源代码；查看修改细节，运行以下代码：

```bash
git diff cb92b24 cec4770
```

---

- `plots/`: 绘图示例
  - `corr_heatmap.ipynb`: 相关性热图绘制
  - `matplotlib_basic.ipynb`: 各种绘图示例

---

- `deprecated/`: 弃用脚本
  - `xsd2vasp.py`: xsd POSCAR 构型文件格式互相转换
