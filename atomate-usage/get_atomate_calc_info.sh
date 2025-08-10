#!/bin/bash

# 批量查看 atomate VASP 计算目录的计算信息
parent_root_dir=$1

# 若 parent... 为空，则 root_dir 为当前目录
if [ -z "$parent_root_dir" ]; then
    # 还有 atomate 计算在运行
    root_dir_list=$(squeue -u yangsl -t RUNNING --format "%.6i %.7T %.8j %.9N %.4C %.11M  %Z" | grep 'FW_job' | awk '{print $7}' | grep "yangsl")
else
    # 所有 atomate 计算已完成
    root_dir_list=$(ls -d $parent_root_dir/*)
fi


for root_dir in $root_dir_list; do
    echo $root_dir | sed "s:/home/yangsl:~:g"

    if [ -z "$parent_root_dir" ]; then
        get_atomate_calc_info.py $root_dir
    else
        get_atomate_calc_info.py $root_dir --flag
    fi

    echo -e "\n"
    printf "%`tput cols`s" | tr ' ' '#'
    echo -e "\n"
done
