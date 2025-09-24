#!/bin/bash

# 获取 AIMD 过程中温度和能量数据并绘制演化图

set -eu

natoms=$(sed -n '7p' POSCAR | awk '{ for(i=1; i<=NF; i++) a+=$i; print a}')
# 获取平均原子能量
grep 'free  energy' OUTCAR | awk -v n=${natoms} '{print $5/n}' > energy.dat
# 获取温度
grep 'T=' OSZICAR | awk '{print $3}' > temperature.dat
# 步数
nsteps=$(grep 'T=' OSZICAR | wc -l)


if hostname | grep -q sjtu; then
  config_path="~/yangsl/scripts/cms-scripts/plots"
else
  config_path="~/scripts/cms-scripts/plots"
fi

# 绘制演化图
cat >> .plot.gnu << EOF
set loadpath "${config_path}"
load "config.gnu"

set output "temperature_energy_evolution.png"

# 横向排列
# set terminal pngcairo size 2000,800
# set multiplot layout 1,2
# 纵向排列
set terminal pngcairo size 1200,1600
set multiplot layout 2,1


# --------- Left Plot: Temperature Evolution ---------
set title "Temperature Evolution"
set xlabel "Ion Step"
set ylabel "Temperature (K)"
set xrange [1:${nsteps}]
plot "temperature.dat" using 0:1 w lp ps 1 pt 3 title "Temperature"

# --------- Right Plot: Energy Evolution ---------
set title "Energy Evolution"
set xlabel "Ion Step"
set ylabel "Energy (eV/atom)"
set xrange [1:${nsteps}]
plot "energy.dat" using 0:1 w lp ps 1 pt 3 title "Energy"

unset multiplot
unset output
EOF

gnuplot .plot.gnu

rm .plot.gnu

echo -e "\nTemperature and energy evolution plot of AIMD saved to temperature_energy_evolution.png."
