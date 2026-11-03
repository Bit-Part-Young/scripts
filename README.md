# CMS Scripts

计算材料领域（Computational Materials Science, CMS）相关实用 Python、Shell 脚本（脚本总数：接近 200）。

内容包括：

- ASE、pymatgen、ovito、phonopy、atomsk、vaspkit 等包/工具使用
- 晶体、缺陷结构建模与操作
- VASP 计算预处理、后处理
- Molecular Dynamics 计算后处理
- 机器学习势函数 MLIPs：NEP、MTP

代码维护：Shenlu Yang, yangsl2126@sjtu.edu.cn。

注意事项：

- 本仓库代码主要依赖 ASE、pymatgen、atomate、PyXtal、spglib、bsym 等 Python 包，以及 atomsk、vaspkit 等工具

- 使用脚本前，请先阅读源码及其中的注释！！！

- 脚本功能主要贴合个人需求及使用，不一定保证通用性，仅供参考！！！

- **具体脚本及功能说明见各子目录的 `README.md`；部分清单可能未及时更新，请以目录中的源码为准**

- 部分脚本的代码参考或拷贝自其他作者，在脚本开头均标明出处；如涉及侵权，请联系作者删除

- 大部分脚本可解析命令行参数，使用 `python xxx.py -h` 或 `./xxx.py -h` 查看帮助信息

- 对于 ipynb 文件，可将其链接复制至 [nbviewer](https://nbviewer.org/) 在线查看

---

## 脚本内容

### ASE 程序使用

- `ase-usage/`: ASE 程序使用，详见 [ase-usage/README.md](ase-usage/README.md)

---

### pymatgen 程序使用

- `pymatgen-usage/`: pymatgen 程序使用，详见 [pymatgen-usage/README.md](pymatgen-usage/README.md)

---

### VASP 相关脚本

- `VASP-scripts/`: VASP 相关脚本，详见 [VASP-scripts/README.md](VASP-scripts/README.md)

---

### phonopy 声子计算

- `phonopy-usage/`: VASP / LAMMPS 声子计算输入生成与声子谱绘图，详见 [phonopy-usage/README.md](phonopy-usage/README.md)

---

### 结构相关脚本

- `structure-scripts/`: 结构相关脚本，详见 [structure-scripts/README.md](structure-scripts/README.md)

---

### NEP 相关脚本

- `NEP`: NEP 相关脚本，详见 [NEP/README.md](NEP/README.md)

---

### MTP 相关脚本

- `MTP`: MTP 相关脚本，详见 [MTP/README.md](MTP/README.md)

---

### MISC 其他

- `ovito-usage/`: ovito Python 使用，详见 [ovito-usage/README.md](ovito-usage/README.md)

- `atomsk-usage/`: atomsk 程序使用，详见 [atomsk-usage/README.md](atomsk-usage/README.md)

---

- `HPC/`: 超算用脚本，详见 [HPC/README.md](HPC/README.md)

---

- `atomate-usage/`: atomate 程序使用，详见 [atomate-usage/README.md](atomate-usage/README.md)

---

- `elastic/`: 弹性相关脚本，详见 [elastic/README.md](elastic/README.md)

---

- `sqsgen-usage/`: sqsgen 程序使用，详见 [sqsgen-usage/README.md](sqsgen-usage/README.md)

---

- `vaspkit-usage/`: vaspkit 程序使用，详见 [vaspkit-usage/README.md](vaspkit-usage/README.md)

---

- `atomkit-usage/`: atomkit 程序使用，详见 [atomkit-usage/README.md](atomkit-usage/README.md)

---

- `image-process/`: 图片处理，详见 [image-process/README.md](image-process/README.md)

---

- `shell-scripts/`: Shell 脚本，详见 [shell-scripts/README.md](shell-scripts/README.md)

---

- `plots/`: 绘图示例，详见 [plots/README.md](plots/README.md)

---

- `misc/`: 其他脚本，详见 [misc/README.md](misc/README.md)
