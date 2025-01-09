#!/bin/bash

# Diamond (100), (110), (111) 表面模型构建

symbol="Si"
# 晶格常数
a=5.43
dup_x=1
dup_y=1
dup_z=1
# 一侧真空层厚度
vacuum=20
# 两侧真空层厚度
vacuum2=$((2 * vacuum))


# (100) 表面
# 或使用 "[01-1]" "[011]" "[100]" 位向
atomsk --create diamond ${a} ${symbol} \
       orient "[010]" "[001]" "[100]" \
       -duplicate ${dup_x} ${dup_y} ${dup_z} \
       -shift 0 0 ${vacuum} \
       -cell add ${vacuum2} z \
       -fractional -sort species pack vasp
    #    lmp

# mv ${symbol}.lmp ${symbol}_100.lmp
mv POSCAR ${symbol}_100.vasp

# (110) 表面
atomsk --create diamond ${a} ${symbol} \
       orient "[1-10]" "[001]" "[110]" \
       -duplicate ${dup_x} ${dup_y} ${dup_z} \
       -shift 0 0 ${vacuum} \
       -cell add ${vacuum2} z \
       -fractional -sort species pack vasp
    #    lmp

# mv ${symbol}.lmp ${symbol}_110.lmp
mv POSCAR ${symbol}_110.vasp

# (111) 表面
atomsk --create diamond ${a} ${symbol} \
       orient "[11-2]" "[-110]" "[111]" \
       -duplicate ${dup_x} ${dup_y} ${dup_z} \
       -shift 0 0 ${vacuum} \
       -cell add ${vacuum2} z \
       -fractional -sort species pack vasp
    #    lmp

# mv ${symbol}.lmp ${symbol}_111.lmp
mv POSCAR ${symbol}_111.vasp

mkdir -p surfaces-atomsk
# mv ${symbol}_100.lmp ${symbol}_110.lmp ${symbol}_111.lmp surfaces-atomsk
mv ${symbol}_100.vasp ${symbol}_110.vasp ${symbol}_111.vasp surfaces-atomsk
