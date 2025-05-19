#!/bin/bash

# 统计 cfg 文件构型帧数及总原子数

set -e
set -u

count_cfg() {

  fcfg=$1

  nframes=$(grep 'BEGIN' ${fcfg} | wc -l)
  natoms=$(sed -n '/Size/{n;p}' ${fcfg} | awk '{sum += $1} END {print sum}')

  echo "${fcfg}: ${nframes} frames, ${natoms} atoms."

}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: $script_name [cfg_fn...]"

  echo -e "\nCount frames and atoms in cfg files."

  echo -e "\nOptions:"
  echo "  -h, --help      Show this help message and exit"
  echo "  cfg_fn          cfg filename"

  echo -e "\nExamples:"
  echo "    ${script_name}                # 当前目录下的所有 cfg 文件"
  echo "    ${script_name} 1.cfg          # 单个 cfg 文件"
  echo "    ${script_name} 1.cfg 2.cfg    # 多个 cfg 文件"

  exit 0
}


#-------------------------------- 主函数 --------------------------------
if [ $# -eq 0 ]; then

  for fcfg in $(ls *.cfg); do
    count_cfg ${fcfg}
  done

elif [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help
  exit 0

else
  for fcfg in "$@"; do
    count_cfg ${fcfg}
  done

fi
