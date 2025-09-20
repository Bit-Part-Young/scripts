#!/bin/bash

# FCC 结构 (100), (110), (111) 表面模型构建

set -eu

surface_fcc() {
  # 设置默认值 元素符号 晶格常数 超胞尺寸 原子层层数 真空层厚度
  surface_type="${surface_type:-100}"
  symbol="${symbol:-Al}"
  a="${a:-4.041}"
  dup_x="${dup_x:-10}"
  dup_y="${dup_y:-10}"
  dup_z="${dup_z:-20}"
  vacuum="${vacuum:-40.0}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")


  #------------------------- (100) 表面 ---------------------------
  # (100) 表面；单胞 2 个原子，2 个原子层
  # 或使用 orient "[010]" "[001]" "[100]"
  if surface_type == "100"; then
    atomsk --create fcc ${a} ${symbol} \
      orient "[01-1]" "[011]" "[100]" \
      -dup ${dup_x} ${dup_y} ${dup_z} \
      -cell add ${vacuum} z \
      -shift 0 0 ${vacuum_half} \
      data.lmp

    echo -e "\n"
    printf "%`tput cols`s" | tr ' ' '#'
    echo -e "\n"

    mv data.lmp 100.lmp


  #------------------------- (110) 表面 ---------------------------
  # (110) 表面；单胞 2 个原子，2 个原子层
  elif surface_type == "110"; then
    atomsk --create fcc ${a} ${symbol} \
      orient "[1-10]" "[001]" "[110]" \
      -dup ${dup_x} ${dup_y} ${dup_z} \
      -cell add ${vacuum} z \
      -shift 0 0 ${vacuum_half} \
      data.lmp

    echo -e "\n"
    printf "%`tput cols`s" | tr ' ' '#'
    echo -e "\n"

    mv data.lmp 110.lmp


  #------------------------- (111) 表面 ---------------------------
  # (111) 表面；单胞 6 个原子，3 个原子层
  elif surface_type == "111"; then
    atomsk --create fcc ${a} ${symbol} \
      orient "[11-2]" "[-110]" "[111]" \
      -dup ${dup_x} ${dup_y} ${dup_z} \
      -cell add ${vacuum} z \
      -shift 0 0 ${vacuum_half} \
      data.lmp

    echo -e "\n"
    printf "%`tput cols`s" | tr ' ' '#'
    echo -e "\n"

    mv data.lmp 111.lmp
  fi

  structure_folder="fcc-surfaces"
  if [[ ! -d ${structure_folder} ]]; then
    mkdir ${structure_folder}
  fi

  mv *.lmp ${structure_folder}

  echo "Total (100), (110), (111) surface models of FCC ${symbol} generated!"
}

#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-st STR] [-e STR] [-lc FLOAT] [-d INT INT INT] [-vac FLOAT]"

  echo -e "\nGenerate FCC (100), (110), (111) surface models for LAMMPS calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -st STR                    surface type (default: 100)"
  echo "    -e STR                     element symbol (default: Al)"
  echo "    -lc FLOAT                  lattice constant (default: 4.041)"
  echo "    -d INT INT INT             x y z dimension (default: 10 10 20)"
  echo "    -vac FLOAT                 vacuum thickness (default: 40.0)"

  echo -e "\nExamples:"
  echo "    ${script_name}"
  echo "    ${script_name} -st 100 -e Al -lc 4.041 -d 10 10 20 -vac 40.0"
  echo "    ${script_name} -st 110 -e Al -lc 4.041 -d 10 10 20 -vac 40.0"
  echo "    ${script_name} -st 111 -e Al -lc 4.041 -d 10 10 10 -vac 40.0"
}


#-------------------------------- 主函数 --------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -st)
      surface_type="$2"
      shift 2
      ;;
    -e)
      symbol="$2"
      shift 2
      ;;
    -lc)
      a="$2"
      shift 2
      ;;
    -d)
      dup_x="$2"
      dup_y="$3"
      dup_z="$4"
      shift 4
      ;;
    -vac)
      vacuum="$2"
      shift 2
      ;;
    -h | --help)
      get_help
      exit 0
      ;;
    *)
      get_help
      exit 1
      ;;
  esac
done

# 调用函数并传递参数
surface_fcc "$@"
