#!/bin/bash

# 获取 VASP 弛豫过程中的能量、原子受力信息（不考虑原子位置方向固定情况，不计算受力的模长）并绘制演化图

# 离子步数
ion_steps=$(grep 'free  ene' OUTCAR | wc -l)

echo "max_x min_x avg_x max_y min_y avg_y max_z min_z avg_z" > force_info.dat
echo "energy_pa" > energy_info.dat

for i in $(seq 1 ${ion_steps}); do

  # 获取平均原子嫩量
  natoms=$(grep 'NIONS' OUTCAR | tail -n 1 | awk '{print $12}')
  energy_pa=$(grep 'free  energy' OUTCAR | sed -n "${i}p" | awk -v n=${natoms} '{print $5/n}')

  echo "${energy_pa}" >> energy_info.dat

  # 获取原子受力信息
  awk '/TOTAL-FORCE/ {flag++} flag==N {print} /total drift/ {if (flag==N) {exit}}' N=${i} OUTCAR | grep -viE 'total|---' | \
  awk '{
      x=$4; y=$5; z=$6;
      if (NR==1) {
          max_x=min_x=avg_x=x;
          max_y=min_y=avg_y=y;
          max_z=min_z=avg_z=z;
      } else {
          if (x>max_x) max_x=x; if (x<min_x) min_x=x; avg_x+=x;
          if (y>max_y) max_y=y; if (y<min_y) min_y=y; avg_y+=y;
          if (z>max_z) max_z=z; if (z<min_z) min_z=z; avg_z+=z;
      }
      n++
  }
  END {
      ax = avg_x / n; ay = avg_y / n; az = avg_z / n
      if (abs(max_x) < 1e-6) max_x = 0.0; if (abs(min_x) < 1e-6) min_x = 0.0; if (abs(ax) < 1e-6) ax = 0.0;
      if (abs(max_y) < 1e-6) max_y = 0.0; if (abs(min_y) < 1e-6) min_y = 0.0; if (abs(ay) < 1e-6) ay = 0.0;
      if (abs(max_z) < 1e-6) max_z = 0.0; if (abs(min_z) < 1e-6) min_z = 0.0; if (abs(az) < 1e-6) az = 0.0;
      print max_x, min_x, ax, max_y, min_y, ay, max_z, min_z, az >> "force_info.dat"
  }
  function abs(val) {
      return val < 0 ? -val : val
  }'

done


# 绘制演化图
cat >> plot_tmp.gnu << EOF
set loadpath "~/scripts/cms-scripts/plots"
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
plot "force_info.dat" using 0:1 w lp ps 4 title "fx_{max}", \
     "force_info.dat" using 0:4 w lp ps 4 title "fy_{max}", \
     "force_info.dat" using 0:7 w lp ps 4 title "fz_{max}"

# --------- Right Plot: Energy Evolution ---------
set title "Energy Evolution"
set xlabel "Ion Step"
set ylabel "Energy (eV/atom)"
set xrange [1:${ion_steps}]
plot "energy_info.dat" using 0:1 w lp ps 4 title "Energy"

unset multiplot
unset output
EOF

gnuplot plot_tmp.gnu

rm plot_tmp.gnu

echo -e "\nEnergy and force evolution plot saved to energy_force_evolution.png."
