# atomate-usage

atomate 程序使用

## 脚本列表

```bash
static.py                            # Diamond Si 静态计算
relax.py                             # Diamond Si 弛豫计算
static_metal.py                      # BCC Nb 静态计算（修改 INCAR、KPOINTS 参数，以适用于金属体系）
relax_metal.py                       # BCC Nb 弛豫计算（修改 INCAR、KPOINTS 参数，以适用于金属体系）
elastic_metal.py                     # FCC Al 弹性常数计算
get_atomate_calc_info.sh             # 批量查看 atomate VASP 计算目录的计算信息
get_atomate_calc_info.py             # 查看 atomate VASP 计算目录的计算信息（配合 get_atomate_calc_info.sh 使用）
atomate_basic.ipynb                  # atomate 基础使用
atomate_db.ipynb                     # 获取使用 atomate 计算并存储到数据库（Mongodb）中的数据
```
