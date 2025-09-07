#!/bin/bash

# 生成 FCC 多晶模型

element="${element:-Al}"
a="${a:-4.05}"


# 生成 poly.txt 配置文件
cat > poly.txt << EOF
box 100 100 300         # 盒子大小
random 6                # 生成 6 个随机取向 & 位置的晶粒
EOF

# 构建单胞
atomsk --create fcc ${a} ${element} ${element}.lmp

# 生成多晶模型
atomsk --polycrystal ${element}.lmp poly.txt -wrap ${element}_polycrystal.lmp

rm ${element}.lmp
