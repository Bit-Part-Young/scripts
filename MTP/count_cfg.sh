#!/bin/bash

# 统计 cfg 构型文件帧数及总原子数

set -e
set -u

for fcfg in $(ls *cfg); do

  nframes=$(grep 'BEGIN' ${fcfg} | wc -l)
  natoms=$(sed -n '/Size/{n;p}' ${fcfg} | awk '{sum += $1} END {print sum}')

  echo "${fcfg}: ${nframes} frames, ${natoms} atoms."

done
