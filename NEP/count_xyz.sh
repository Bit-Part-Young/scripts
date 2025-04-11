#!/bin/bash

# 统计 extxyz 构型文件帧数及总原子数
# reference: https://github.com/Kick-H/For_gpumd/blob/master/NEP_related/Count/count_xyz.sh

set -e
set -u

for fxyz in $(ls *xyz); do

  nlines=$(wc -l ${fxyz} | cut -d ' ' -f1)
  natoms=$(grep -i pbc ${fxyz} | wc -l)

  echo -e "${fxyz}: \c"
  echo ${nlines} ${natoms} | awk '{print $2 " frames, " $1-$2*2 " atoms."}'

done
