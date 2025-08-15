#!/bin/bash

# 获取当前 目录/子目录 下的 VASP 静态/单点能计算数据


#-------------------------------- 获取 VASP 静态/单点能计算数据 --------------------------------
get_scf_data() {

  # 输出表头
  echo "|-----------------------------------------------------------------------------------------------------------------------"
  printf "%21s %8s %10s %13s %15s %11s %16s %13s\n" "Folder" "natoms" "Ele_Step" "Energy" "dE" "Energy_pa" "TimeCost" "State"
  echo "|-----------------------------------------------------------------------------------------------------------------------"

  unfinished_dirs=()

  for dir in $(echo "$@" | tr ' ' '\n' | sort -n); do
    if [ -d "$dir" ]; then
      oszicar_fn="${dir}/OSZICAR"
      outcar_fn="${dir}/OUTCAR"
      if [ -f "${oszicar_fn}" ]; then

        # 若目录层级数大于 2，则获取最后 3 层目录
        if [[ $(echo "${dir}" | awk -F'/' '{print NF}') -gt 2 ]]; then
          # dir=$(echo "${dir}" | awk -F'/' '{print $(NF-1)"/"$NF}')
          dir=$(echo "${dir}" | awk -F'/' '{print $(NF-2)"/"$(NF-1)"/"$NF}')
        else
          dir=$(basename "${dir}")
        fi

        # 只取后 21 个字符
        dir=$(echo "${dir}" | rev | cut -c -21 | rev)

        if grep -q 'DAV' "${oszicar_fn}"; then
          natoms=$(grep 'NIONS' "${outcar_fn}" | tail -n 1 | awk '{print $12}')
          nsteps=$(grep 'DAV' "${oszicar_fn}" | tail -n 1 | awk '{print $2}')

          time=$(grep 'Total CPU time used' "${outcar_fn}" | awk '{print $6}')
          time=${time%.*}
          if [[ ${time} -gt 0 ]]; then
            state="Completed"

            dE=""

            # 运行结束时需使用 F= 行所对应的 E0 数据
            energy=$(grep 'F=' "${oszicar_fn}" | tail -n 1 | awk '{printf "%.6f", $5}')

          elif [[ ${time} -eq 0 ]]; then
            state="Running"

            time=$(grep 'LOOP:' "${outcar_fn}" | awk '{sum+=$7} END {print sum}')
            time=${time%.*}

            dE=$(grep 'DAV' "${oszicar_fn}" | tail -n 1 | awk '{print $4}')

            energy=$(grep 'DAV' "${oszicar_fn}" | tail -n 1 | awk '{printf "%.6f", $3}')
          fi
          energy_pa=$(awk "BEGIN { print ${energy} / ${natoms} }")

          hours=$((time / 3600))
          minutes=$(((time % 3600) / 60))
          seconds=$((time % 60))
          time_info=$(printf "%02dh %02dm %02ds" "${hours}" "${minutes}" "${seconds}")


          printf "%21s %8s %10s %13s %15s %11s %16s %13s\n" "${dir}" "${natoms}" "${nsteps}" "${energy}" "${dE}" "${energy_pa}" "${time_info}" "${state}"
        fi
      fi
    fi
  done

  echo "|-----------------------------------------------------------------------------------------------------------------------"

}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: $script_name [dirs...]"

  echo -e "\nGet electronic steps, energy and time cost data in VASP scf calculation dirs with table format."

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
      get_scf_data "$item"
    fi
  done

else
  get_scf_data "$@"
fi
