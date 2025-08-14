#!/bin/bash

: '
Generate ISF, ESF/TWIN configurations of FCC {111}<112> slip system for VASP GSFE calculation.

上半部分原子层向右滑移 a/√6，形成 ISF 模型；下半部分原子层（去除第一个原子层）向左滑移 a/√6，形成 ESF/TWIN 模型（实际只为 ESF 模型）。

reference: Materials Today Communications 43 (2025) 111639, https://doi.org/10.1016/j.mtcomm.2025.111639
'


#-------------------------------- Generate stacking fault configurations --------------------------------
sf_generation() {
  # 设置默认值
  element="${element:-Al}"
  a="${a:-4.041}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  dup_z="${dup_z:-4}"
  num_interval="${num_interval:-10}"
  vacuum="${vacuum:-15.0}"


  interlayer_distance=$(echo "scale=15; ${a}*sqrt(3)/3" | bc -l)
  length_z=$(echo "scale=15; ${a}*sqrt(3)*${dup_z}" | bc -l)

  length_intrinsic=$(echo "scale=15; 0.985*(${length_z}/2)" | bc -l)
  length_twin=$(echo "scale=15; 0.985*(${length_z}/2-${interlayer_distance})" | bc -l)

  # calculate the slip step length
  b_interval=$(echo "scale=15; ${a}/sqrt(6)/${num_interval}" | bc -l)


  #------------------------- 初始位向超胞构建 ---------------------------
  printf "%`tput cols`s" | tr ' ' '#'
  echo -e "\n\nInitial configuration:\n"

  atomsk --create fcc ${a} ${element} \
    orient "[1-10]" "[11-2]" "[111]" \
    -duplicate ${dup_x} ${dup_y} ${dup_z} \
    -fix x -fix y \
    -sort species pack \
    -fractional -ow vasp

  echo


  #-------------------------     创建存储结构文件夹     ---------------------------
  folders=("ISF" "TWIN")

  for folder in "${folders[@]}"; do
      if [ -d "${folder}" ]; then
          rm -rf "${folder}"
          mkdir "${folder}"
      else
          mkdir "${folder}"
      fi
  done


  #------------------------- ISF 模型 ---------------------------
  cd ISF

  for i in $(seq 0 1 ${num_interval}); do
    disp=$(echo "scale=15; ${b_interval}*${i}" | bc)

    atomsk ../POSCAR \
      -shift above ${length_intrinsic} z 0.0 ${disp} 0.0 \
      -cell add ${vacuum} z \
      -fractional -wrap -ow vasp

    mv POSCAR "POSCAR.${i}"

    echo

  done

  printf "%`tput cols`s" | tr ' ' '#'
  echo

  cd ..


  #------------------------- TWIN 模型 ---------------------------
  cd TWIN

  for i in $(seq 0 1 ${num_interval}); do
    disp=$(echo "scale=15; -1*${b_interval}*${i}" | bc)

    atomsk ../ISF/POSCAR.${num_interval} \
      -shift below ${length_twin} z 0 ${disp} 0 \
      -fractional -wrap -ow vasp

    mv POSCAR "POSCAR.${i}"

    echo

  done

  printf "%`tput cols`s" | tr ' ' '#'
  echo

  cd ..
}


#-------------------------------- Get help --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-e symbol] [-lc lattice_constant] [-d dup_x dup_y dup_z] [-vac vacuum] [-ni num_interval]"

  echo -e "\nGenerate ISF, ESF/TWIN configurations of FCC {111}<112> slip system for VASP GSFE calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -e element                 element symbol (default: Al)"
  echo "    -lc lattice_constant       lattice constant (default: 4.041)"
  echo "    -d dup_x dup_y dup_z       duplicate system in the three directions (default: 1 1 6)"
  echo "    -vac vacuum                vacuum thickness (default: 15.0)"
  echo "    -ni num_interval           number of intervals (default: 10)"

  echo -e "\nExamples:"
  echo "    Default settings: ${script_name}"
  echo "    For pure system: ${script_name} -e Al -lc 4.041 -d 1 1 4 -vac 15.0 -ni 10"
}


#-------------------------------- Main function --------------------------------
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
    -ni)
      num_interval="$2"
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


# call the function and pass parameters
sf_generation "$@"
