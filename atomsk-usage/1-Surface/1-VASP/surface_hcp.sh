#!/bin/bash

# HCP (001), (100) 表面模型构建

symbol="Ti"
# 晶格常数
a=2.928
c=4.640
dup_x=1
dup_y=1
dup_z=6
# 一侧真空层厚度
vacuum=5
# 两侧真空层厚度
vacuum2=$((2 * vacuum))


# (001) 或 Basal 表面
atomsk --create hcp ${a} ${c} ${symbol} \
       orient "[11-20]" "[-1100]" "[0001]" \
       -duplicate ${dup_x} ${dup_y} ${dup_z} \
       -shift 0 0 ${vacuum} \
       -cell add ${vacuum2} z \
       -fractional -sort species pack vasp

mv POSCAR ${symbol}_001.vasp

# (100) 或 Prismatic 表面
atomsk --create hcp ${a} ${c} ${symbol} \
       orient "[11-20]" "[0001]" "[-1100]" \
       -duplicate ${dup_x} ${dup_y} 3 \
       -shift 0 0 ${vacuum} \
       -cell add ${vacuum2} z \
       -fractional -sort species pack vasp

mv POSCAR ${symbol}_100.vasp

structure_folder="0-Structures"
mkdir -p ${structure_folder}
mv ${symbol}_001.vasp ${symbol}_100.vasp ${structure_folder}
