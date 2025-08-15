#!/bin/bash

# HCP Basal / {0001}, Prismatic / {11-20} , Pyramidal I / {10-11}, Pyramidal II / {11-22} 表面模型构建


surface_hcp() {
  # 设置默认值 元素符号 晶格常数 超胞尺寸 原子层层数 真空层厚度
  symbol="${symbol:-Ti}"
  a="${a:-2.928}"
  c="${c:-4.640}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  num_layer="${num_layer:-12}"
  vacuum="${vacuum:-15.0}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")


  #------------------------- Basal / (0001) 表面 ---------------------------
  # 默认扩胞 1x1x6
  dup_z=$(awk "BEGIN {print ${num_layer}/2}")
  atomsk --create hcp ${a} ${c} ${symbol} \
    orient "[11-20]" "[-1100]" "[0001]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
    -sort species pack -fractional vasp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv POSCAR POSCAR.basal


  #------------------------- Prismatic / {11-20} 表面 ---------------------------
  # 默认扩胞 1x1x3
  dup_z=$(awk "BEGIN {print ${num_layer}/4}")
  atomsk --create hcp ${a} ${c} ${symbol} \
    orient "[11-20]" "[0001]" "[-1100]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
    -sort species pack -fractional vasp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv POSCAR POSCAR.prismatic


  #------------------------- Pyramidal I / {10-11} 表面 ---------------------------
  # 默认扩胞 1x1x2
  dup_z=$(awk "BEGIN {print ${num_layer}/6}")
  atomsk --create hcp ${a} ${c} ${symbol} \
    orient "[-2110]" "[1-213]" "[01-11]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
    -sort species pack -fractional vasp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv POSCAR POSCAR.pyramidalI


  #------------------------- Pyramidal II / {11-22} 表面 ---------------------------
  # 默认扩胞 1x1x1
  dup_z=$(awk "BEGIN {print ${num_layer}/12}")
  atomsk --create hcp ${a} ${c} ${symbol} \
    orient "[-1-123]" "[-1100]" "[11-22]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
    -sort species pack -fractional vasp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv POSCAR POSCAR.pyramidalII


  #------------------------- 创建存储结构文件的目录 ---------------------------
  structure_folder="hcp-surfaces"
  if [[ ! -d ${structure_folder} ]]; then
    mkdir ${structure_folder}
  fi

  mv POSCAR.basal POSCAR.prismatic POSCAR.pyramidalI POSCAR.pyramidalII ${structure_folder}

  echo "Total HCP Basal / {0001}, Prismatic / {11-20}, Pyramidal I / {10-11}, Pyramidal II / {11-22} surface models of HCP ${symbol} generated!"
}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-e symbol] [-lc lattice_constant] [-d dup_x dup_y num_layer] [-vac vacuum]"

  echo -e "\nGenerate HCP Basal / {0001}, Prismatic / {11-20}, Pyramidal I / {10-11}, Pyramidal II / {11-22} surface models for VASP calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -e element                 element symbol (default: Ti)"
  echo "    -lc lattice_constant       lattice constant (default: 2.928 4.640)"
  echo "    -d dup_x dup_y num_layer   x y dimension, layers of z direction (default: 1 1 12)"
  echo "    -vac vacuum                vacuum thickness (default: 15.0)"

  echo -e "\nExamples:"
  echo "    ${script_name}"
  echo "    ${script_name} -e Ti -lc 2.928 4.640 -d 1 1 12 -vac 15.0"
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
      c="$3"
      shift 3
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
surface_hcp "$@"
