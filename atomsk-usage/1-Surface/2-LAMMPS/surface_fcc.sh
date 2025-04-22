#!/bin/bash

# FCC 结构 (100), (110), (111) 表面模型构建

# [ ] 是否要固定中间原子层？


surface_fcc() {
  # 设置默认值 元素符号 晶格常数 超胞尺寸 原子层层数 真空层厚度
  symbol="${symbol:-Al}"
  a="${a:-4.041}"
  dup_x="${dup_x:-5}"
  dup_y="${dup_y:-5}"
  dup_z="${dup_z:-10}"
  vacuum="${vacuum:-40}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")


  # (100) 表面；单胞 2 个原子，2 个原子层
  # 或使用 orient "[010]" "[001]" "[100]"
  atomsk --create fcc ${a} ${symbol} \
    orient "[01-1]" "[011]" "[100]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -fix x -fix y \
    -cell add ${vacuum} z \
    -shift 0 0 ${vacuum_half} \
    data.lmp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv data.lmp 100.lmp


  # (110) 表面；单胞 2 个原子，2 个原子层
  atomsk --create fcc ${a} ${symbol} \
    orient "[1-10]" "[001]" "[110]" \
    -fix x -fix y \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -cell add ${vacuum} z \
    -shift 0 0 ${vacuum_half} \
    data.lmp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv data.lmp 110.lmp


  # (111) 表面；单胞 6 个原子，3 个原子层
  atomsk --create fcc ${a} ${symbol} \
    orient "[11-2]" "[-110]" "[111]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -fix x -fix y \
    -cell add ${vacuum} z \
    -shift 0 0 ${vacuum_half} \
    data.lmp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv data.lmp 111.lmp

  structure_folder="fcc-surfaces"
  if [[ ! -d ${structure_folder} ]]; then
    mkdir ${structure_folder}
  fi

  mv 100.lmp 110.lmp 111.lmp ${structure_folder}

  echo "Total (100), (110), (111) surface models of FCC ${symbol} generated!"
}

#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-e symbol] [-lc lattice_constant] [-d dup_x dup_y dup_z] [-vac vacuum]"

  echo -e "\nGenerate FCC (100), (110), (111) surface models for LAMMPS calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -e element                 element symbol (default: Al)"
  echo "    -lc lattice_constant       lattice constant (default: 4.041)"
  echo "    -d dup_x dup_y dup_z       x y dimension, layers of z direction (default: 5 5 10)"
  echo "    -vac vacuum                vacuum thickness (default: 40)"

  echo -e "\nExamples:"
  echo "    ${script_name}"
  echo "    ${script_name} -e Al -lc 4.041 -d 5 5 10 -vac 40"
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
