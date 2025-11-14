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

echo -e "\nCurrent nsteps: ${nsteps}.\n"
echo "Statistic(last half of data):"

# 统计 energy.dat 和 temperature.dat 两个文件后一半行数的数据的平均值和标准差
awk '
{ x[NR]=$1 }                                  # 读入所有行
END {
    start = int(NR/2) + 1                      # 后半部分起始行
    n = NR - start + 1

    # 平均值
    for(i=start;i<=NR;i++){ sum += x[i] }
    mean = sum / n

    # 标准差
    for(i=start;i<=NR;i++){ ss += (x[i]-mean)^2 }
    std = sqrt(ss / n)

    mean = sprintf("%.2f", mean)
    std = sprintf("%.3f", std)
    print "Energy mean:", mean, "std:", std
}' energy.dat

awk '
{ x[NR]=$1 }
END {
    start = int(NR/2) + 1
    n = NR - start + 1

    for(i=start;i<=NR;i++){ sum += x[i] }
    mean = sum / n

    for(i=start;i<=NR;i++){ ss += (x[i]-mean)^2 }
    std = sqrt(ss / n)

    mean = sprintf("%.1f", mean)
    std = sprintf("%.1f", std)
    print "Temperature mean:", mean, "std:", std
}' temperature.dat


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
