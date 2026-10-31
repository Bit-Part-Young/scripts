# phonopy-usage

phonopy 声子谱计算（VASP 有限位移 / DFPT，LAMMPS 势函数校验）的输入生成与绘图脚本，需要 phonopy ≥ 4。各脚本用 `-h` 查看参数。

| 文件 | 功能 |
|---|---|
| `phonon_vasp.sh` | 生成 VASP 声子计算输入：换成原胞、扩胞，`-m fd` 每个位移一个 `disp-XXX/`，`-m dfpt` 一个 `DFPT/`；写好 INCAR、POTCAR 和作业脚本 |
| `disp_supercell_generation.py` | 用 phonopy API 生成位移超胞，`--interface vasp / lammps`；LAMMPS 格式用于势函数声子谱校验 |
| `band_{bcc,fcc,hcp}.conf` | phonopy 命令行算能带用的 k 路径（Γ-H-P-Γ-N、Γ-X-U\|K-Γ-L、Γ-K-M-Γ-A），与 `plot_phonon_bs.py` 的预设一致 |
| `plot_phonon_bs.py` | 主要绘图脚本：读含力的 phonopy yaml，按预设 k 路径计算并绘制声子谱；可叠加多组数据、实验散点，用 YAML 配置画多子图 |
| `plot_band_yaml.py` | 直接画 phonopy 命令行输出的 `band.yaml`，可叠加多个文件；画图风格同 `plot_phonon_bs.py` |
| `plot_mp_phonon.py` | 从 Materials Project 下载声子能带（dfpt / pheasy）并作图，需 `PMG_MAPI_KEY` |

## 典型流程

VASP 有限位移法：

```bash
phonon_vasp.sh -d 4 4 4                    # 生成 disp-XXX/，提交计算
phonopy-init -f disp-*/vasprun.xml         # 收集力，得到 FORCE_SETS
phonopy --config band_bcc.conf -s          # 得到 band.yaml
plot_band_yaml.py band.yaml -o Nb.png      # 作图
```

DFPT：`phonon_vasp.sh -m dfpt -d 3 3 3`，计算完成后 `phonopy-init --fc DFPT/vasprun.xml` 得到 FORCE_CONSTANTS，
conf 中改用 `FORCE_CONSTANTS = READ`。

LAMMPS 势函数：`disp_supercell_generation.py POSCAR -d 6 6 4 --interface lammps`，
对每个 `supercell-XXX` 做单点计算并 dump 力（格式见脚本说明），再 `phonopy-init --lammps -f lammps_forces_*.0`，其余同上。

需要和其他势函数、实验值对比时，用 `plot_phonon_bs.py`，见其文件头的说明。

依赖：`phonon_vasp.sh` 调用 cms-scripts 中的 `get_psp2.py`、`check_vaspi`（`VASP-scripts/`）和 `slurm_generation.py`（`HPC/`）。
