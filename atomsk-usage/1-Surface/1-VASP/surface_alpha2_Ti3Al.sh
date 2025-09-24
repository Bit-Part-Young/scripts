#!/bin/bash

set -eu

# alpha2-Ti3Al Basal / {0001} 表面模型构建


surface_alpha2_Ti3Al() {
  # 设置默认值 晶格常数 超胞尺寸 原子层层数 真空层厚度
  a="${a:-5.753}"
  c="${c:-4.650}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  num_layer="${num_layer:-4}"
  vacuum="${vacuum:-15.0}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")

  initial_structure_fn="alpha2_Ti3Al.vasp"
  bx=$(echo "scale=10; -1*${a}/2" | bc)
  by=$(echo "scale=10; ${a}*sqrt(3)/2" | bc)

  cat > ${initial_structure_fn} << EOF
# alpha2-Ti3Al unit cell
1.0
   ${a}    0.0000000000    0.0000000000
   ${bx}   ${by}           0.0000000000
   0.0000000000    0.0000000000    ${c}
Ti Al
6 2
direct
  0.8303216930  0.6606433860  0.2500000000 Ti
  0.3393566140  0.1696783070  0.2500000000 Ti
  0.8303216930  0.1696783070  0.2500000000 Ti
  0.1696783070  0.3393566140  0.7500000000 Ti
  0.6606433860  0.8303216930  0.7500000000 Ti
  0.1696783070  0.8303216930  0.7500000000 Ti
  0.3333333333  0.6666666667  0.2500000000 Al
  0.6666666667  0.3333333333  0.7500000000 Al
EOF


  #------------------------- Basal / (0001) 表面 ---------------------------
  atomsk ${initial_structure_fn} \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -cell add ${vacuum} z -shift 0 0 ${vacuum_half} \
    -sort species pack -fractional vasp

  rm ${initial_structure_fn}

  echo -e "\n"
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n"

  mv POSCAR POSCAR.basal

  #------------------------- 创建存储结构文件的目录 ---------------------------
  structure_folder="alpha2-Ti3Al-surfaces"
  if [[ ! -d ${structure_folder} ]]; then
    mkdir ${structure_folder}
  fi

  mv POSCAR.basal ${structure_folder}

  echo "alpha2-Ti3Al Basal / {0001} surface model generated!"
}



#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-lc lattice_constant] [-d dup_x dup_y num_layer] [-vac vacuum]"

  echo -e "\nGenerate alpha2-Ti3Al Basal / {0001} surface model for VASP calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -lc lattice_constant       lattice constant (default: 5.753 4.650)"
  echo "    -d dup_x dup_y num_layer   x y dimension, layers of z direction (default: 1 1 4)"
  echo "    -vac vacuum                vacuum thickness (default: 15.0)"

  echo -e "\nExamples:"
  echo "    ${script_name}"
  echo "    ${script_name} -lc 5.753 4.650 -d 1 1 4 -vac 15.0"
}


#-------------------------------- 主函数 --------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
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
surface_alpha2_Ti3Al "$@"
