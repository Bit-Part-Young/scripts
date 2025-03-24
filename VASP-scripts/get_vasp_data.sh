#!/bin/bash


#-------------------------------- 获取 VASP 数据 --------------------------------
get_vasp_data() {

  # 输出表头
  echo "|-----------------------------------------------------------------------------------"
  echo "|        Folder         |   Step   |       Energy       |  Energy_pa  |  TimeCost  |"
  echo "|-----------------------------------------------------------------------------------"

  unfinished_dirs=()

  for dir in "$@"; do
    if [ -d "$dir" ]; then
      oszicar_fn="${dir}/OSZICAR"
      outcar_fn="${dir}/OUTCAR"
      if [ -f "$oszicar_fn" ]; then
        dir=$(basename "${dir}")
        if grep -q 'F=' "$oszicar_fn"; then
          energy=$(grep 'F=' "$oszicar_fn" | tail -n 1 | awk '{print $5}')
          natoms=$(grep 'NIONS' "$outcar_fn" | tail -1 | awk '{print $12}')
          energy_pa=$(awk "BEGIN { print ${energy} / ${natoms} }")

          ionstep=$(grep 'F=' "$oszicar_fn" | tail -n 1 | awk '{print $1}')

          time=$(grep 'Total CPU time used' "$outcar_fn" | awk '{print $6}')

          printf "| %-21s | %-8s | %-18s | %-11s | %-10s |\n" "${dir}" "${ionstep}" "${energy}" "${energy_pa}" "${time}"
          echo "|-----------------------------------------------------------------------------------"
        else
          unfinished_dirs+=("${dir}")
        fi
      fi
    fi
  done

  # 统一输出未完成第一个离子步的 VASP 计算目录
  if [[ ${#unfinished_dirs[@]} -gt 0 ]]; then
    echo -e "\nVASP calculation dirs did not finish the first ion step:"
    for dir in "${unfinished_dirs[@]}"; do
      echo "  $dir"
    done
  fi
}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: $script_name [options] [dirs...]"

  echo -e "\nGet ion steps, energy and time cost data in VASP calculation dirs with table format."

  echo -e "\nOptions:"
  echo "  -h, --help      Show this help message and exit"
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
      get_vasp_data "$item"
    fi
  done

else
  get_vasp_data "$@"
fi
