# atomsk-usage

atomsk 程序使用

## 脚本列表

```bash
# 1-Surface/ 目录
surface_bcc.sh                       # BCC (100), (110), (111) 表面模型构建
surface_fcc.sh                       # FCC (100), (110), (111) 表面模型构建
surface_hcp.sh                       # HCP Basal, Prismatic , Pyramidal I, Pyramidal II 表面模型构建
surface_diamond.sh                   # Diamond (100), (110), (111) 表面模型构建
surface_alpha2_Ti3Al.sh              # alpha2-Ti3Al Basal / {0001} 表面模型构建
surface_gamma_TiAl.sh                # gamma-TiAl (111) 表面模型构建；重要近似条件 c≈a

# 2-Stacking-Fault/ 目录
sf_fcc.sh                            # 生成 FCC {111}<112> 滑移系的 ISF, ESF 和 TWIN 层错构型
sf_fcc2.sh                           # 生成 FCC {111}<112> 滑移系的 ISF, ESF/TWIN 层错构型（效果与 sf_fcc.sh 的前两步一致）
twin_fcc.sh                          # FCC 孪晶模型构建

1-Al-polycrystal/                    # FCC Al 多晶模型构建
2-Fe-Cr-Ni-polycrystal/              # Fe-Cr-Ni 多晶模型构建
```
