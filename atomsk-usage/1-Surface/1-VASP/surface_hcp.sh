#!/bin/bash

# HCP Basal / (0001), Prismatic / (-1100), Pyramidal I / (01-11), Pyramidal II / (11-22) 表面模型构建

set -eu

surface_hcp() {
  # 设置默认值 元素符号 晶格常数 超胞尺寸 原子层层数 真空层厚度
  surface_type="${surface_type:-basal}"
  symbol="${symbol:-Ti}"
  a="${a:-2.928}"
  c="${c:-4.640}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  dup_z="${dup_z:-6}"
  vacuum="${vacuum:-15.0}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")


  #------------------------- Basal / (0001) 表面 ---------------------------
  # 默认扩胞 1x1x6
  if [[ ${surface_type} == "basal" ]]; then
    atomsk --create hcp ${a} ${c} ${symbol} \
      orient "[11-20]" "[-1100]" "[0001]" \
      -dup ${dup_x} ${dup_y} ${dup_z} \
      -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
      -sort species pack -fractional vasp

    echo -e "\n"
    printf "%`tput cols`s" | tr ' ' '#'
    echo -e "\n"

    mv POSCAR POSCAR.basal


  #------------------------- Prismatic / (-1100) 表面 ---------------------------
  # 默认扩胞 1x1x3
  elif [[ ${surface_type} == "prismatic" ]]; then
    atomsk --create hcp ${a} ${c} ${symbol} \
      orient "[11-20]" "[0001]" "[-1100]" \
      -dup ${dup_x} ${dup_y} ${dup_z} \
      -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
      -sort species pack -fractional vasp

    echo -e "\n"
    printf "%`tput cols`s" | tr ' ' '#'
    echo -e "\n"

    mv POSCAR POSCAR.prismatic


  #------------------------- Pyramidal I / (01-11) 表面 ---------------------------
  # 默认扩胞 1x1x2
  elif [[ ${surface_type} == "pyramidalI" ]]; then
    atomsk --create hcp ${a} ${c} ${symbol} \
      orient "[-2110]" "[1-213]" "[01-11]" \
      -dup ${dup_x} ${dup_y} ${dup_z} \
      -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
      -sort species pack -fractional vasp

    echo -e "\n"
    printf "%`tput cols`s" | tr ' ' '#'
    echo -e "\n"

    mv POSCAR POSCAR.pyramidalI


  #------------------------- Pyramidal II / (11-22) 表面 ---------------------------
  # 默认扩胞 1x1x1
  elif [[ ${surface_type} == "pyramidalII" ]]; then
    atomsk --create hcp ${a} ${c} ${symbol} \
      orient "[-1-123]" "[-1100]" "[11-22]" \
      -dup ${dup_x} ${dup_y} ${dup_z} \
      -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
      -sort species pack -fractional vasp

    echo -e "\n"
    printf "%`tput cols`s" | tr ' ' '#'
    echo -e "\n"

    mv POSCAR POSCAR.pyramidalII
  fi

  #------------------------- 创建存储结构文件的目录 ---------------------------
  structure_folder="hcp-surfaces"
  if [[ ! -d ${structure_folder} ]]; then
    mkdir ${structure_folder}
  fi

  mv POSCAR.* ${structure_folder}

  echo "Total HCP Basal / (0001), Prismatic / (-1100), Pyramidal I / (01-11), Pyramidal II / (11-22) surface models of HCP ${symbol} generated!"
}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-st STR] [-e STR] [-lc FLOAT FLOAT] [-d INT INT INT] [-vac FLOAT]"

  echo -e "\nGenerate HCP Basal / (0001), Prismatic / (-1100), Pyramidal I / (01-11), Pyramidal II / (11-22) surface models for VASP calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -st STR                    surface type (default: basal)"
  echo "    -e STR                     element symbol (default: Ti)"
  echo "    -lc FLOAT FLOAT            lattice constant (default: 2.928 4.640)"
  echo "    -d INT INT INT             x y z dimension (default: 1 1 6)"
  echo "    -vac FLOAT                 vacuum thickness (default: 15.0)"

  echo -e "\nExamples:"
  echo "    ${script_name}"
  echo "    ${script_name} -st basal -e Ti -lc 2.928 4.640 -d 1 1 6 -vac 15.0"
  echo "    ${script_name} -st prismatic -e Ti -lc 2.928 4.640 -d 1 1 3 -vac 15.0"
  echo "    ${script_name} -st pyramidalI -e Ti -lc 2.928 4.640 -d 1 1 2 -vac 15.0"
  echo "    ${script_name} -st pyramidalII -e Ti -lc 2.928 4.640 -d 1 1 1 -vac 15.0"
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
      c="$3"
      shift 3
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
surface_hcp "$@"
