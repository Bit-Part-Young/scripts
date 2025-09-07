#!/bin/bash

# 生成 Fe-Cr-Ni 多晶模型


# 生成 poly.txt 配置文件
cat > poly.txt << EOF
box 100 100 300
random 12
EOF

# 构建 Ni 单胞
atomsk --create fcc 3.48 Ni -duplicate 5 5 5 Ni.xsf

# 置换 Ni 为 Fe
atomsk Ni.xsf -select random 30% Ni -substitute Ni Fe Fe_Ni.xsf

# 置换 Ni 为 Cr
atomsk Fe_Ni.xsf -select random 20% Ni -substitute Ni Cr Fe_Cr_Ni.xsf

# 生成多晶模型
atomsk --polycrystal Fe_Cr_Ni.xsf poly.txt -wrap incoloy_poly.lmp

rm *.xsf
