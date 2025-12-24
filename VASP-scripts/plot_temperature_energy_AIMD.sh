#!/bin/bash

# 获取 AIMD 过程中温度和能量数据并绘制演化图

set -eu


#----------------------- 统计文件后一半行数的数据的平均值和标准差 ---------------------------
# $1=数据文件名, $2=标签（Temperature / Energy）, $3=小数点输出格式
calc_stats() {
    awk -v label="$2" -v fmt="$3" '
    { x[NR]=$1 }
    END {
        start = int(NR/2) + 1
        n = NR - start + 1

        for(i=start;i<=NR;i++){ sum += x[i] }
        mean = sum / n

        for(i=start;i<=NR;i++){ ss += (x[i]-mean)^2 }
        std = sqrt(ss / n)

        mean = sprintf(fmt, mean)
        std = sprintf(fmt, std)
        print label " mean:", mean, "std:", std
    }' "$1"
}


#----------------------- 获取温度、能量、压强随步数演化的数据 ---------------------------
natoms=$(sed -n '7p' POSCAR | awk '{ for(i=1; i<=NF; i++) a+=$i; print a}')
# 获取平均原子能量
grep 'free  energy' OUTCAR | awk -v n=${natoms} '{print $5/n}' > energy.dat
# 获取温度
grep 'T=' OSZICAR | awk '{print $3}' > temperature.dat
# 获取压强
grep 'external pressure' OUTCAR | awk '{print $4}' > pressure.dat
# 步数
nsteps=$(grep 'T=' OSZICAR | wc -l)

echo -e "\nCurrent nsteps: ${nsteps}.\n"
echo "Statistic(last half of data):"
calc_stats energy.dat "Energy" "%.2f"
calc_stats temperature.dat "Temperature" "%.1f"
calc_stats pressure.dat "Pressure" "%.1f"


# ----------------------- 绘制演化图 -------------------------
if hostname | grep -q sjtu; then
  config_path="~/yangsl/scripts/cms-scripts/plots"
else
  config_path="~/scripts/cms-scripts/plots"
fi

cat > .plot.gnu << EOF
set loadpath "${config_path}"
load "config.gnu"

set output "temperature_energy_evolution.png"

# 横向排列
# set terminal pngcairo size 2000,800
# set multiplot layout 1,2
# 纵向排列
set terminal pngcairo size 1200,1600
set multiplot layout 2,1


# --------- upper Plot: Temperature Evolution ---------
set title "Temperature Evolution"
set xlabel "Ion Step"
set ylabel "Temperature (K)"
set xrange [1:${nsteps}]
plot "temperature.dat" using 0:1 w lp ps 1 pt 3 title "Temperature"

# --------- lower Plot: Energy Evolution ---------
set title "Energy Evolution"
set xlabel "Ion Step"
set ylabel "Energy (eV/atom)"
set xrange [1:${nsteps}]
plot "energy.dat" using 0:1 w lp ps 1 pt 3 title "Energy"

unset multiplot
unset output
unset terminal
unset title


# --------- Pressure Evolution Plot ---------
set terminal pngcairo size 1000,800
set output "pressure_evolution.png"

set xlabel "Ion Step"
set ylabel "Pressure (GPa)"
set xrange [1:${nsteps}]
plot "pressure.dat" using 0:(\$1/10) w lp ps 1 pt 3 title "Pressure"

EOF

gnuplot .plot.gnu

rm .plot.gnu

echo -e "\nTemperature, energy, and pressure evolution plot of AIMD saved to temperature_energy_evolution.png and pressure_evolution.png."
