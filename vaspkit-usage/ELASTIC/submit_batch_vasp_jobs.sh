#!/bin/bash

# 批量提交由 vaspkit 生成的 能量-应变法 弹性常数计算任务

flag=$1

root_path=$(pwd)
for cij in $(ls -F | grep /$ | grep ^C); do
  cd ${root_path}/$cij
  for s in strain_*; do
    cd ${root_path}/$cij/$s

    # 删除 KPOINTS 符号链接
    if [[ -L "KPOINTS" ]]; then
      rm KPOINTS
    fi

    cp ../../job.slurm .


    if [[ ${flag} == "check" ]]; then
      check_vaspi
    elif [[ ${flag} == "submit" ]]; then
      sbatch job.slurm
    fi

  done
done
