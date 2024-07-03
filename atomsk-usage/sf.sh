#!/bin/bash

# 构建层错

# 构建含 111 滑移系 {111}<110> 的 Al bulk
atomsk --create fcc 4.046 Al orient [-110] [111] [11-2] -duplicate 1 8 1 Al_cell.xsf
# 沿滑移方向移动 1 埃
atomsk Al_cell.xsf -shift above 0.5*box Y 1.0 0.0 0.0 -wrap Al_sf.xsf

# 格式转换
atomsk Al_cell.xsf -sort species pack -fractional vasp
mv POSCAR Al_cell.vasp
atomsk Al_SF.xsf -sort species pack -fractional vasp
mv POSCAR Al_SF.vasp

rm *.xsf
