#!/bin/bash

: '
FCC 晶体的 ISF（本征堆垛层错）, ESF（非本征堆垛层错） 和 TWIN（孪晶）层错模型构建

reference:
- https://mp.weixin.qq.com/s/o0ldM-87tGulOEIBFXSsSw
- https://mp.weixin.qq.com/s/5Mk37UB0SwwbLAw7FaZRZg
'


#-------------------------     创建存储结构文件夹     ---------------------------
folders=("ISF" "ESF" "TWIN")

for folder in "${folders[@]}"; do
    if [ -d "${folder}" ]; then
        rm -rf "${folder}"
        mkdir "${folder}"
    else
        mkdir "${folder}"
    fi
done


#-------------------------     定义变量     ---------------------------
# 晶格常数
lat=3.53
# 超胞
dup_x=1
dup_y=2
dup_z=10
# 真空层厚度
vacuum=20
# 两侧真空层厚度
vacuum2=$((2 * vacuum))

# z 方向的原子层间距
gap_z=$(awk "BEGIN {print ${lat}*sqrt(3)/3}")
# z 方向长度
len_z=$(awk "BEGIN {print ${lat}*sqrt(3)*${dup_z}+${vacuum2}}")

# z 方向长度一半
len_i=$(awk "BEGIN {print ${len_z}/2}")
# z 方向长度一半 + 1 * z 方向原子层间距
len_e=$(awk "BEGIN {print ${len_i}+${gap_z}}")
# z 方向长度一半 + 2 * z 方向原子层间距
len_t=$(awk "BEGIN {print ${len_e}+${gap_z}}")

# a*[112]/6 分 40 次滑移
ndisp=40
# 沿 [11-2] 方向每次移动的距离
disp=$(awk "BEGIN {print ${lat}/sqrt(6)/${ndisp}}")


#-------------------------     初始位向超胞构建     ---------------------------
atomsk --create fcc ${lat} Ni \
       orient "[1-10]" "[11-2]" "[111]" \
       -duplicate ${dup_x} ${dup_y} ${dup_z} \
       -shift 0 0 ${vacuum} \
       -cell add ${vacuum2} z \
       -wrap Ni_init.lmp


#-------------------------     ISF 模型     ---------------------------
atomsk Ni_init.lmp Ni_0_isf.lmp
for i in {1..40..1}; do
  atomsk Ni_$((i-1))_isf.lmp \
         -shift above ${len_i} z 0 ${disp} 0 \
         -wrap Ni_${i}_isf.lmp
done


#-------------------------     ESF 模型     ---------------------------
atomsk Ni_40_isf.lmp Ni_0_esf.lmp
for i in {1..40..1}; do
  atomsk Ni_$((i-1))_esf.lmp \
         -shift above ${len_e} z 0 ${disp} 0 \
         -wrap Ni_${i}_esf.lmp
done


#-------------------------     TWIN 模型     ---------------------------
atomsk Ni_40_esf.lmp Ni_0_twin.lmp
for i in {1..40..1}; do
  atomsk Ni_$((i-1))_twin.lmp \
         -shift above ${len_t} z 0 ${disp} 0 \
         -wrap Ni_${i}_twin.lmp
done


# 移动构型文件至对应的文件夹
mv Ni_*_isf.lmp  ISF
mv Ni_*_esf.lmp  ESF
mv Ni_*_twin.lmp TWIN