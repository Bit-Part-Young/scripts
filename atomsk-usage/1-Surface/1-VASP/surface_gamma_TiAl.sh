#!/bin/bash

# gamma-TiAl (111) 表面模型构建；重要近似条件 c≈a


surface_gamma_TiAl() {
  # 设置默认值 晶格常数 超胞尺寸 原子层层数 真空层厚度
  a="${a:-3.993}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  dup_z="${dup_z:-4}"
  vacuum="${vacuum:-15.0}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")


  #------------------------- (111) 表面 ---------------------------
  # (111) 表面；单胞 6 个原子，3 个原子层
  atomsk --create fcc ${a} Al Ti \
    orient "[-110]" "[11-2]" "[111]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
    -sort species pack -fractional vasp

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv POSCAR POSCAR.111

  structure_folder="gamma-TiAl-surfaces"
  if [[ ! -d ${structure_folder} ]]; then
    mkdir ${structure_folder}
  fi

  mv POSCAR.111 ${structure_folder}

  echo "gamma-TiAl (111) surface model generated!"
}

#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-lc lattice_constant] [-d dup_x dup_y dup_z] [-vac vacuum]"

  echo -e "\nGenerate gamma-TiAl (111) surface model for VASP calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -lc lattice_constant       lattice constant (default: 3.993)"
  echo "    -d dup_x dup_y dup_z       x y dimension, layers of z direction (default: 1 1 4)"
  echo "    -vac vacuum                vacuum thickness (default: 15.0)"

  echo -e "\nExamples:"
  echo "    ${script_name}"
  echo "    ${script_name} -lc 3.993 -d 1 1 4 -vac 15.0"
}


#-------------------------------- 主函数 --------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
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
surface_gamma_TiAl "$@"
