#!/bin/bash

# 生成 FCC Al 多晶模型

# 构建 Al 单胞
atomsk --create fcc 4.05 Al Al.lmp

# 多晶模型生成
atomsk --polycrystal Al.lmp poly.txt -wrap Al_polycrystal.lmp

rm Al.lmp
