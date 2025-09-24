#!/bin/bash

set -eu

# Diamond (100), (110), (111) 表面模型构建


surface_diamond() {
  # 设置默认值 元素符号 晶格常数 超胞尺寸 原子层层数 真空层厚度
  symbol="${symbol:-Si}"
  a="${a:-5.43}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  # [ ] 待确定 num_layer 和 dup_z
  num_layer="${num_layer:-12}"
  vacuum="${vacuum:-15}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")


  # (100) 表面
  # 或使用 "[01-1]" "[011]" "[100]" 位向
  dup_z=$(awk "BEGIN {print ${num_layer}/2}")
  atomsk --create diamond ${a} ${symbol} \
    orient "[010]" "[001]" "[100]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -shift 0 0 ${vacuum} \
    -cell add ${vacuum_half} z \
    -sort species pack \
    -fractional vasp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv POSCAR POSCAR.100

  # (110) 表面
  atomsk --create diamond ${a} ${symbol} \
    orient "[1-10]" "[001]" "[110]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -shift 0 0 ${vacuum} \
    -cell add ${vacuum_half} z \
    -sort species pack \
    -fractional vasp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv POSCAR POSCAR.110

  # (111) 表面
  atomsk --create diamond ${a} ${symbol} \
    orient "[11-2]" "[-110]" "[111]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -shift 0 0 ${vacuum} \
    -cell add ${vacuum_half} z \
    -sort species pack \
    -fractional vasp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv POSCAR POSCAR.111

  structure_folder="diamond-surfaces"
  if [[ ! -d ${structure_folder} ]]; then
    mkdir ${structure_folder}
  fi

  mv POSCAR.100 POSCAR.110 POSCAR.111 ${structure_folder}

  echo "Total (100), (110), (111) surface models of Diamond ${symbol} generated!"

}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-e symbol] [-lc lattice_constant] [-d dup_x dup_y num_layer] [-vac vacuum]"

  echo -e "\nGenerate Diamond (100), (110), (111) surface models for VASP calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -e element                 element symbol (default: Si)"
  echo "    -lc lattice_constant       lattice constant (default: 5.43)"
  echo "    -d dup_x dup_y num_layer   x y dimension, layers of z direction (default: 1 1 12)"
  echo "    -vac vacuum                vacuum thickness (default: 15)"

  echo -e "\nExamples:"
  echo "    ${script_name}"
  echo "    ${script_name} -e Si -lc 5.43 -d 1 1 12 -vac 15"
}


#-------------------------------- 主函数 --------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
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
      num_layer="$4"
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
surface_diamond "$@"
