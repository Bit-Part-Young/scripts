#!/bin/bash

# 输出多个 VASP 计算目录下的输出数据

# 输出表头
echo "|-----------------------------------------------------------"
echo "|        Folder         |      Energy       |   TimeCost   |"
echo "|-----------------------------------------------------------"

for item in $(ls); do
  if [ -d "$item" ]; then
    if [ -f "${item}/OSZICAR" ]; then
        energy=$(grep 'F=' "${item}/OSZICAR" | tail -n 1 | awk '{print $5}')
        time=$(grep 'Total CPU time used' "${item}/OUTCAR" | awk '{print $6}')

        printf "| %-21s | %-17s | %-12s |\n" "${item}" "${energy}" "${time}"
        echo "|-----------------------------------------------------------"
    fi
  fi
done
