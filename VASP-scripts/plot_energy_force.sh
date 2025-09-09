#!/bin/bash

# 获取 VASP 弛豫过程中的能量、原子受力（考虑原子位置方向固定情况，计算受力的模长）信息并绘制演化图

set -e

awk '{if($4=="F"||$4=="T") print $4,$5,$6}' POSCAR > temp.fixed

# 离子步数
ion_steps=$(grep 'free  ene' OUTCAR | wc -l)

echo "force_max" > force_info.dat
echo "energy_pa" > energy_info.dat

for i in $(seq 1 ${ion_steps}); do

  # 获取平均原子嫩量
  natoms=$(grep 'NIONS' OUTCAR | tail -n 1 | awk '{print $12}')
  energy_pa=$(grep 'free  energy' OUTCAR | sed -n "${i}p" | awk -v n=${natoms} '{print $5/n}')

  echo "${energy_pa}" >> energy_info.dat

  # 获取原子受力信息
  awk '/TOTAL-FORCE/ {flag++} flag==N {print} /total drift/ {if (flag==N) {exit}}' N=${i} OUTCAR \
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
      awk ' BEGIN { max = 0 } FNR==NR { flag[FNR] = $0; next }
      {
          split(flag[FNR], f, /[ \t]+/)
          for (i = 1; i <= NF; i++) {
              if (f[i] == "F") $i = 0
          }
          sum = 0
          for (i = 1; i <= NF; i++) sum += $i * $i
          norm = sqrt(sum)
          if (norm > max) max = norm
      }
      END { print max }' temp.fixed temp.force >> force_info.dat
  fi

done

rm temp.fixed temp.force


if hostname | grep -q sjtu; then
  config_path="~/yangsl/scripts/cms-scripts/plots"
else
  config_path="~/scripts/cms-scripts/plots"
fi

# 绘制演化图
cat >> .plot.gnu << EOF
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
fi

rm .plot.gnu

echo -e "\nEnergy and force evolution plot saved to energy_force_evolution.png."
