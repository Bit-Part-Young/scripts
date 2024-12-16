#!/bin/bash

# 以 TUI 表格形式输出当前多个 VASP 计算目录下的离子步、能量、耗时数据

# 输出表头
echo "|-----------------------------------------------------------------------------------"
echo "|        Folder         |   Step   |       Energy       |  Energy_pa  |  TimeCost  |"
echo "|-----------------------------------------------------------------------------------"

for item in $(ls); do
  if [ -d "$item" ]; then
    if [ -f "${item}/OSZICAR" ]; then
        energy=$(grep 'F=' "${item}/OSZICAR" | tail -n 1 | awk '{print $5}')
        natoms=$(grep 'NIONS' "${item}/OUTCAR" | tail -1 | awk '{print $12}')
        energy_pa=$(awk "BEGIN { print ${energy} / ${natoms} }")

        ionstep=$(grep 'F=' "${item}/OSZICAR" | tail -n 1 | awk '{print $1}')

        time=$(grep 'Total CPU time used' "${item}/OUTCAR" | awk '{print $6}')

        printf "| %-21s | %-8s | %-18s | %-11s | %-10s |\n" "${item}" "${ionstep}" "${energy}" "${energy_pa}" "${time}"
        echo "|-----------------------------------------------------------------------------------"
    fi
  fi
done
