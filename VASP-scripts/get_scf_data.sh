#!/bin/bash

# 获取当前 目录/子目录 下的 VASP 单点能计算数据


#-------------------------------- 获取 VASP 单点能计算数据 --------------------------------
get_scf_data() {

  # 输出表头
  echo "|------------------------------------------------------------------------------------------------------------"
  printf "%21s %8s %17s %13s %15s %11s   %-13s\n" "Folder" "natoms" "Electronic_Step" "Energy" "dE" "Energy_pa" "TimeCost"
  echo "|------------------------------------------------------------------------------------------------------------"

  unfinished_dirs=()

  for dir in $(echo "$@" | tr ' ' '\n' | sort -n); do
    if [ -d "$dir" ]; then
      oszicar_fn="${dir}/OSZICAR"
      outcar_fn="${dir}/OUTCAR"
      if [ -f "$oszicar_fn" ]; then
        dir=$(basename "${dir}")

        if grep -q 'DAV' "$oszicar_fn"; then
          electronic_step=$(grep 'DAV' "$oszicar_fn" | tail -n 1 | awk '{print $2}')
          energy=$(grep 'DAV' "$oszicar_fn" | tail -n 1 | awk '{printf "%.6f", $3}')
          natoms=$(grep 'NIONS' "$outcar_fn" | tail -1 | awk '{print $12}')
          energy_pa=$(awk "BEGIN { print ${energy} / ${natoms} }")

          time=$(grep 'Total CPU time used' "$outcar_fn" | awk '{print $6}')
          time=${time%.*}
          if [[ ${time} -gt 0 ]]; then
            hours=$((time / 3600))
            minutes=$(((time % 3600) / 60))
            seconds=$((time % 60))
            time=$(printf "%02dh %02dm %02ds" "$hours" "$minutes" "$seconds")

            dE=""
          elif [[ ${time} -eq 0 ]]; then
            time="Still Running"

            dE=$(grep 'DAV' "$oszicar_fn" | tail -n 1 | awk '{print $4}')
          fi

          printf "%21s %8s %17s %13s %15s %11s   %-13s\n" "${dir}" "${natoms}" "${electronic_step}" "${energy}" "${dE}" "${energy_pa}" "${time}"
        fi
      fi
    fi
  done

  echo "|------------------------------------------------------------------------------------------------------------"

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
