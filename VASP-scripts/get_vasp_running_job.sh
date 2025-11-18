#!/bin/bash

set -eu

# 获取正在运行的 VASP 作业信息（scf 和 opt 类型）

cal_type=$1

platform=$(hostname)

if [[ ${platform} =~ "sjtu" ]]; then
    dir_list=($(squeue -t RUNNING --format "%.6i %.7T %.8j %.9N %.4C %.11M  %Z" | grep 'VASP' | awk '{print $7}' | grep "yangsl"))
else
    dir_list=($(squeue -u yangsl -t RUNNING --format "%.6i %.7T %.8j %.9N %.4C %.11M  %Z" | grep 'VASP' | awk '{print $7}' | grep "yangsl"))
fi

if [ "$cal_type" == "scf" ]; then
    get_scf_data.sh ${dir_list[@]}
elif [ "$cal_type" == "opt" ]; then
    get_opt_data.sh ${dir_list[@]}
fi
