#!/bin/bash

# 统计 xyz 文件中的 element 字段类别及其数量
# reference: https://github.com/Kick-H/For_gpumd/blob/master/NEP_related/Count/count_xyz.sh

set -e
set -u

count_xyz_type() {

  fxyz=$1

  echo -e "\n$fxyz count info:"

  count_total=0
  elements=$(grep -o 'element=[^ ]*' ${fxyz} | awk -F= '{print $2}' | sort | uniq)
  for element in ${elements[@]}; do
    echo -e "${element}: \c"
    count=$(grep -E element=${element} ${fxyz} | wc -l)
    echo ${count}
    count_total=$((count_total + count))
  done
  echo "Total: ${count_total}"

}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: $script_name [xyz_fn...]"

  echo -e "\nCount element types in xyz files."

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

  for fxyz in $(ls *.xyz); do
    count_xyz_type ${fxyz}
  done

elif [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help
  exit 0

else
  for fxyz in "$@"; do
    count_xyz_type ${fxyz}
  done

fi
