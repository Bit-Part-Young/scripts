#!/bin/bash

set -eu

# 统计 xyz 文件构型帧数及总原子数
# reference: https://github.com/Kick-H/For_gpumd/blob/master/NEP_related/Count/count_xyz.sh


count_xyz() {

  xyz_fn=$1

  if [ ! -f ${xyz_fn} ]; then
    echo "Error: ${xyz_fn} not found!"
    exit 1
  fi

  nlines=$(wc -l ${xyz_fn} | cut -d ' ' -f1)
  natoms=$(grep -i pbc ${xyz_fn} | wc -l)

  echo -e "${xyz_fn}: \c"
  echo ${nlines} ${natoms} | awk '{print $2 " frames, " $1-$2*2 " atoms."}'

}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: $script_name [xyz_fn...]"

  echo -e "\nCount frames and atoms in xyz files."

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
    count_xyz ${xyz_fn}
  done

elif [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help
  exit 0

else
  for xyz_fn in "$@"; do
    count_xyz ${xyz_fn}
  done

fi
