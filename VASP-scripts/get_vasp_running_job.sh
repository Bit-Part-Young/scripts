#!/bin/bash

# 获取正在运行的 VASP 作业信息

type=$1

root_dir_list=($(squeue -u yangsl -t RUNNING --format "%.6i %.7T %.8j %.9N %.4C %.11M  %Z" | grep 'VASP' | awk '{print $7}' | grep "yangsl"))

if [ "$type" == "scf" ]; then
    get_scf_data.sh ${root_dir_list[@]}
elif [ "$type" == "opt" ]; then
    get_opt_data.sh ${root_dir_list[@]}
fi
