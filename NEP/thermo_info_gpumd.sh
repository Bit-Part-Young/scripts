#!/usr/bin/bash

# 从 GPUMD 输出的 thermo.out 文件中格式化输出热力学信息

#-------------------------------- 从 thermo.out 文件中格式化输出信息 --------------------------------
get_thermo_info() {
  info_type="${1:-all}"

  # 跳过前 18 行
  tmp_file="thermo.out.tmp"
  awk 'NR>18' thermo.out > ${tmp_file}

  column=$(awk 'END { print NF }' ${tmp_file})

  if [[ ${info_type} == "all" ]]; then
    if [[ ${column} -eq 12 ]]; then
      printf "%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n" "T" "Ek" "Ep" "Px" "Py" "Pz" "Pyz" "Pxz" "Pxy" "Lx" "Ly" "Lz"
    elif [[ ${column} -eq 18 ]]; then
      printf "%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n" "T" "Ek" "Ep" "Px" "Py" "Pz" "Pyz" "Pxz" "Pxy" "ax" "ay" "az" "bx" "by" "bz" "cx" "cy" "cz"
    fi

    awk '{for(i=1; i<=NF; i++) printf "%-15.6f ", $i; print ""}' ${tmp_file}
  elif [[ ${info_type} == "temp" ]]; then
    printf "%-15s %-15s %-15s\n" "T" "Ek" "Ep"

    awk '{for(i=1; i<=3; i++) printf "%-15.6f ", $i; print ""}' ${tmp_file}
  elif [[ ${info_type} == "press" ]]; then
    printf "%-15s %-15s %-15s %-15s %-15s %-15s\n" "Px" "Py" "Pz" "Pyz" "Pxz" "Pxy"

    awk '{for(i=4; i<=9; i++) printf "%-15.6f ", $i; print ""}' ${tmp_file}
  elif [[ ${info_type} == "lattice" ]]; then
    if [[ ${column} -eq 12 ]]; then
      printf "%-15s %-15s %-15s\n" "Lx" "Ly" "Lz"

      awk '{for(i=NF-2; i<=NF; i++) printf "%-15.6f ", $i; print ""}' ${tmp_file}
    elif [[ ${column} -eq 18 ]]; then
      printf "%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n" "ax" "ay" "az" "bx" "by" "bz" "cx" "cy" "cz"

      awk '{for(i=NF-8; i<=NF; i++) printf "%-15.6f ", $i; print ""}' ${tmp_file}
    fi
  fi

  rm -f ${tmp_file}
}

#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [info_type]"

  echo -e "\nPrint thermo info from gpumd thermo.out file."

  echo -e "\nOptions:"
  echo "    -h, --help    show this help message and exit"
  echo "    info_type     info type, temp, press, lattice"

  echo -e "\nExamples:"
  echo "    ${script_name}           # 输出所有信息"
  echo "    ${script_name} all       # 输出所有信息"
  echo "    ${script_name} temp      # 输出温度、动能、势能"
}

#-------------------------------- 主函数 --------------------------------
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help
  exit 0
elif [[ $# -eq 0 ]]; then
  get_thermo_info
elif [[ $# -eq 1 ]]; then
  get_thermo_info "$@"
fi