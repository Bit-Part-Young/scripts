#!/bin/bash

# 获取当前 目录/子目录 下的 VASP 弛豫计算数据


#-------------------------------- 获取 VASP 弛豫计算数据 --------------------------------
get_opt_data() {

  # 输出表头
  echo "|---------------------------------------------------------------------------------------------"
  printf "%21s %10s %13s %11s %16s %13s\n" "Folder" "Ion_Step" "Energy" "Energy_pa" "TimeCost" "State"
  echo "|---------------------------------------------------------------------------------------------"

  unfinished_dirs=()

  for dir in $(echo "$@" | tr ' ' '\n' | sort -n); do
    if [ -d "$dir" ]; then
      oszicar_fn="${dir}/OSZICAR"
      outcar_fn="${dir}/OUTCAR"
      incar_fn="${dir}/INCAR"
      if [ -f "$oszicar_fn" ]; then

        # 若目录层级数大于 2，则获取最后 3 层目录
        if [[ $(echo "${dir}" | awk -F'/' '{print NF}') -gt 2 ]]; then
          # dir=$(echo "${dir}" | awk -F'/' '{print $(NF-1)"/"$NF}')
          dir=$(echo "${dir}" | awk -F'/' '{print $(NF-2)"/"$(NF-1)"/"$NF}')
        else
          dir=$(basename "${dir}")
        fi

        # 只取后 21 个字符
        dir=$(echo "${dir}" | rev | cut -c -21 | rev)

        if grep -q 'F=' "${oszicar_fn}"; then

          # 能量、温度数据获取获取
          if grep -q -E '^TEBEG' "${incar_fn}"; then
            energy=$(grep 'F=' "${oszicar_fn}" | tail -n 1 | awk '{printf "%.6f", $9}')
            # temperature=$(grep 'TEBEG' "${outcar_fn}" | awk -F';' '{print $1}' | awk '{print $3}')
          else
            energy=$(grep 'F=' "${oszicar_fn}" | tail -n 1 | awk '{printf "%.6f", $5}')
          fi

          natoms=$(grep 'NIONS' "${outcar_fn}" | tail -n 1 | awk '{print $12}')
          energy_pa=$(awk "BEGIN { print ${energy} / ${natoms} }")

          nsteps=$(grep 'F=' "${oszicar_fn}" | tail -n 1 | awk '{print $1}')

          # 获取耗时数据
          time=$(grep 'Total CPU time used' "${outcar_fn}" | awk '{print $6}')
          time=${time%.*}
          if [[ ${time} -gt 0 ]]; then
            state="Completed"
          elif [[ ${time} -eq 0 ]]; then
            time=$(grep 'LOOP+' "${outcar_fn}" | awk '{sum+=$7} END {print sum}')
            time=${time%.*}
            state="Running"
          fi

          hours=$((time / 3600))
          minutes=$(((time % 3600) / 60))
          seconds=$((time % 60))
          time_info=$(printf "%02dh %02dm %02ds" "${hours}" "${minutes}" "${seconds}")

          printf "%21s %10s %13s %11s %16s %13s\n" "${dir}" "${nsteps}" "${energy}" "${energy_pa}" "${time_info}" "${state}"

        else
          unfinished_dirs+=("${dir}")
        fi
      fi
    fi
  done

  echo "|---------------------------------------------------------------------------------------------"

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

  echo -e "\nUsage: $script_name [dirs...]"

  echo -e "\nGet ion steps, energy and time cost data in VASP relaxation calculation dirs with table format."

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
      get_opt_data "$item"
    fi
  done

else
  get_opt_data "$@"
fi
