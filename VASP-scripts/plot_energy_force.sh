#!/bin/bash

# 获取 VASP 弛豫过程中的能量、原子受力（考虑原子位置方向固定情况，计算受力的模长）信息并绘制演化图

set -eu

# 只读取 Selective dynamics 行之后的坐标行，避免第 6 行元素行（如含 F 元素）被误识别
awk 'NR >= 7 && /^[[:space:]]*[Ss]/ {sd = 1; next} sd && toupper($4) ~ /^[TF]$/ {print $4,$5,$6}' POSCAR > temp.fixed

# 离子步数
ion_steps=$(grep 'free  ene' OUTCAR | wc -l)

echo "force_max" > force_info.dat
echo "energy_pa" > energy_info.dat

for i in $(seq 1 ${ion_steps}); do

  # 获取平均原子嫩量
  natoms=$(grep 'NIONS' OUTCAR | tail -n 1 | awk '{print $12}')
  energy_pa=$(grep 'free  energy' OUTCAR | sed -n "${i}p" | awk -v n=${natoms} '{print $5/n}')

  echo "${energy_pa}" >> energy_info.dat

  # 获取原子受力信息；同时将该离子步的晶格（前 3 列正格矢、后 3 列倒格矢）写入 temp.latt
  awk '/direct lattice vectors/ {getline l1; getline l2; getline l3}
       /TOTAL-FORCE/ {flag++; if (flag==N) print l1 "\n" l2 "\n" l3 > "temp.latt"}
       flag==N {print} /total drift/ {if (flag==N) {exit}}' N=${i} OUTCAR \
  | grep -viE 'total|---' | awk '{print $4,$5,$6}' > temp.force

 if [[ ! -s "temp.fixed" ]]; then
      awk 'BEGIN { max = 0 } {
          sum = 0
          for (i = 1; i <= NF; i++) sum += $i * $i
          norm = sqrt(sum)
          if (norm > max) max = norm
      }
      END { print max }' temp.force >> force_info.dat
  else
      # 与 VASP 源码 dyna.F 中 SET_SELECTED_FORCES_ZERO 一致：
      # Selective dynamics 的 T/F 对应晶格矢量 a、b、c 方向（而非笛卡尔 x、y、z 轴），
      # 先用倒格矢 B 将受力转为分量 c_j = F·b_j（KARDIR），将固定方向分量置 0，
      # 再用正格矢 A 转回笛卡尔坐标 F' = sum_j c_j·a_j（DIRKAR）；
      # A、B 取自当前离子步（temp.latt），正交晶胞下等价于直接将对应的笛卡尔分量置 0
      awk ' BEGIN { max = 0 }
      FILENAME == ARGV[1] {
          for (k = 1; k <= 3; k++) { A[FNR, k] = $k; B[FNR, k] = $(k + 3) }
          next
      }
      FILENAME == ARGV[2] { flag[FNR] = $0; next }
      {
          split(flag[FNR], f, /[ \t]+/)
          for (j = 1; j <= 3; j++) {
              c[j] = 0
              for (k = 1; k <= 3; k++) c[j] += $k * B[j, k]
              if (toupper(f[j]) == "F") c[j] = 0
          }
          sum = 0
          for (k = 1; k <= 3; k++) {
              g = 0
              for (j = 1; j <= 3; j++) g += c[j] * A[j, k]
              sum += g * g
          }
          norm = sqrt(sum)
          if (norm > max) max = norm
      }
      END { print max }' temp.latt temp.fixed temp.force >> force_info.dat
  fi

done

rm -f temp.fixed temp.force temp.latt


if hostname | grep -q sjtu; then
  config_path="~/yangsl/scripts/cms-scripts/plots"
else
  config_path="~/scripts/cms-scripts/plots"
fi

# 绘制演化图
# 使用 > 覆盖，避免上次异常退出残留的 .plot.gnu 被重复追加
cat > .plot.gnu << EOF
set loadpath "${config_path}"
load "config.gnu"

set output "energy_force_evolution.png"

# 横向排列
set terminal pngcairo size 2000,800
set multiplot layout 1,2
# 纵向排列
# set terminal pngcairo size 1200,1600
# set multiplot layout 2,1


# --------- Left Plot: Force Evolution ---------
set title "Force Evolution"
set xlabel "Ion Step"
set ylabel "Force (eV/Å)"
set xrange [1:${ion_steps}]
plot "force_info.dat" using 0:1 w lp ps 3 pt 3 title "f_{max}"

# --------- Right Plot: Energy Evolution ---------
set title "Energy Evolution"
set xlabel "Ion Step"
set ylabel "Energy (eV/atom)"
set xrange [1:${ion_steps}]
plot "energy_info.dat" using 0:1 w lp ps 3 pt 3 title "Energy"

unset multiplot
unset output
EOF

if [[ "${ion_steps}" -gt 1 ]]; then
  gnuplot .plot.gnu
  echo -e "\nEnergy and force evolution plot saved to energy_force_evolution.png."
else
  echo -e "\nOnly ${ion_steps} ion step, skip plotting (data in force_info.dat, energy_info.dat)."
fi

rm .plot.gnu
