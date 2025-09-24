#!/bin/bash

set -eu

# 统计 VASP 计算目录的离子步步数、计算耗时、能量信息

outcar_path=$1

if [ -z "$outcar_path" ]; then
    outcar_path="."
fi

outcar_fn="${outcar_path}/OUTCAR"
nsteps=$(grep 'LOOP+' $outcar_fn | wc -l)
nsteps=$(echo "${nsteps}" | sed 's/ //g')

energy=$(grep '  energy  without entropy=' "${outcar_fn}" | awk '{printf "%.5f\n", $7}' | tail -n 1)

time_cost=$(grep 'LOOP+' "${outcar_fn}" | awk '{sum+=$7} END {print sum}')

hour=$(echo "${time_cost} / 3600" | bc)
minute=$(echo "${time_cost} % 3600 / 60" | bc)
second=$(echo "${time_cost} % 60" | bc | awk '{printf "%.0f\n", $1}')

echo -e "\n${nsteps} inonic steps, current energy: ${energy} eV, time cost: ${hour}h ${minute}m ${second}s"
