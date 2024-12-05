#!/bin/bash

# FCC (100), (110), (111) 表面模型构建

symbol="Cu"
# 晶格常数
a=3.615
# 一侧真空层厚度
vacuum=20
# 两侧真空层厚度
vacuum2=$((2 * vacuum))


# (100) 表面
atomsk --create fcc ${a} ${symbol} orient "[010]" "[001]" "[100]" -duplicate 10 10 6 -shift 0 0 ${vacuum} -cell add ${vacuum2} z lmp

mv ${symbol}.lmp ${symbol}_100.lmp

# (110) 表面
atomsk --create fcc ${a} ${symbol} orient "[1-10]" "[001]" "[110]" -duplicate 10 10 6 -shift 0 0 ${vacuum} -cell add ${vacuum2} z lmp

mv ${symbol}.lmp ${symbol}_110.lmp

# (111) 表面
atomsk --create fcc ${a} ${symbol} orient "[11-2]" "[-110]" "[111]" -duplicate 10 10 6 -shift 0 0 ${vacuum} -cell add ${vacuum2} z lmp

mv ${symbol}.lmp ${symbol}_111.lmp

mkdir -p surfaces-atomsk
mv ${symbol}_100.lmp ${symbol}_110.lmp ${symbol}_111.lmp surfaces-atomsk
