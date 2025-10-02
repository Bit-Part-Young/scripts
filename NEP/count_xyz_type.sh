#!/bin/bash

set -eu

# 统计 xyz 文件中的 element & group 字段类别及其数量
# reference: https://github.com/Kick-H/For_gpumd/blob/master/NEP_related/Count/count_xyz.sh


count_xyz_type() {

  xyz_fn=$1

  if [ ! -f ${xyz_fn} ]; then
    echo "Error: ${xyz_fn} not found!"
    exit 1
  fi

  echo -e "\n${xyz_fn} info:"

  # ----------------------------- 统计 element 类型 -----------------------------
  echo -e "\nElement type count:"
  count_total=0
  elements=$(grep -o 'element=[^ ]*' ${xyz_fn} | awk -F= '{print $2}' | sort | uniq)
  for element in ${elements[@]}; do
    echo -e "${element}: \c"
    count=$(grep -E element=${element} ${xyz_fn} | wc -l)
    echo ${count}
    count_total=$((count_total + count))
  done
  echo "Total: ${count_total}."

  # ----------------------------- 统计 group 类型 -----------------------------
  echo -e "\nGroup type count:"
  count_total=0
  groups=$(grep -o 'group=[^ ]*' ${xyz_fn} | awk -F= '{print $2}' | sort | uniq)
  for group in ${groups[@]}; do
    echo -e "${group}: \c"
    count=$(grep -E group=${group} ${xyz_fn} | wc -l)
    echo ${count}
    count_total=$((count_total + count))
  done
  echo "Total: ${count_total}."
}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: $script_name [xyz_fn...]"

  echo -e "\nCount element & group type in extxyz file."

  echo -e "\nOptions:"
  echo "  -h, --help      Show this help message and exit"
  echo "  xyz_fn          xyz filename"

  echo -e "\nExamples:"
  echo "    ${script_name}                # 当前目录下的所有 xyz 文件"
  echo "    ${script_name} 1.xyz          # 单个 xyz 文件"
  echo "    ${script_name} 1.xyz 2.xyz    # 多个 xyz 文件"

  exit 0
}


#-------------------------------- 主函数 --------------------------------
if [ $# -eq 0 ]; then

  for xyz_fn in $(ls *.xyz); do
    count_xyz_type ${xyz_fn}
  done

elif [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help
  exit 0

else
  for xyz_fn in "$@"; do
    count_xyz_type ${xyz_fn}
  done

fi
