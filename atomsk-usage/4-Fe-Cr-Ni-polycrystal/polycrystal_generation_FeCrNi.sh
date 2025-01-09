#!/bin/bash

# 生成 Fe-Cr-Ni 多晶模型
# -properties props.txt *.lmp 无效？
# 转而使用 xyz 文件构型，可不用 -properties 选项

# 构建 Ni 单胞
atomsk --create fcc 3.48 Ni -duplicate 5 5 5 Ni.xyz

# 置换 Ni 为 Fe
atomsk Ni.xyz -select random 30% Ni -substitute Ni Fe Fe_Ni.xyz

# 置换 Ni 为 Cr
atomsk Fe_Ni.xyz -select random 20% Ni -substitute Ni Cr Fe_Cr_Ni.xyz

# 多晶模型生成
atomsk --polycrystal Fe_Cr_Ni.xyz poly.txt -wrap incoloy_poly.lmp

rm Ni.xyz Fe_Ni.xyz
