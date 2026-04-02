#!/bin/bash

set -eu

# 检查当前 目录/子目录 下的 VASP 是否有计算错误


#-------------------------------- 检查 VASP 是否有计算错误 --------------------------------
check_vasp_errors() {

  errors_type=(
    "Error EDDDAV"
    "WARNING in EDDRMM:"
    "WARNING: Sub-Space-Matrix"
    "BRMIX: very serious problems"
    "LAPACK: Routine ZPOTRF failed"
  )

  echo
  for error_type in "${errors_type[@]}"; do
    echo -n "\"${error_type}\" check: "
    for dir in $(echo "$@" | tr ' ' '\n' | sort -n); do
      if [[ -d "$dir" ]]; then
           grep -q "${error_type}" "${dir}"/*.out 2>/dev/null && echo -n "${dir} "
      fi
    done

    echo -e "\n"
  done

}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: $script_name [dirs...]"

  echo -e "\nCheck errors in VASP calculation dirs."

  echo -e "\nOptions:"
  echo "  -h, --help      show this help message and exit"
  echo "  dirs            VASP calculation dirs"

  echo -e "\nExamples:"
  echo "    ${script_name} \${PWD}       # 当前目录"
  echo "    ${script_name} test1        # 单个目录"
  echo "    ${script_name} test1 test2  # 多个目录"
  echo "    ${script_name} .            # 当前目录下的所有子目录"

  exit 0
}


#-------------------------------- 主函数 --------------------------------
if [ $# -eq 0 ]; then
  get_help
  exit 0

elif [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help

elif [[ "$1" == "." ]]; then
  for item in ./*; do
    if [ -d "$item" ]; then
      check_vasp_errors "$item"
    fi
  done

else
  check_vasp_errors "$@"
fi
