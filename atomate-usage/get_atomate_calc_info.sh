#!/bin/bash

# 批量查看 atomate VASP 计算目录的计算信息

root_dir_list=$(squeue -u yangsl -t RUNNING --format "%.6i %.7T %.8j %.9N %.4C %.11M  %Z" | grep 'FW_job' | awk '{print $6}' | grep "yangsl")

for root_dir in $root_dir_list; do
    echo $root_dir | sed "s:/home/yangsl:~:g"

    get_atomate_calc_info.py $root_dir

    echo -e "\n"
    printf "%`tput cols`s" | tr ' ' '#'
    echo -e "\n"
done