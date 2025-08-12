#!/bin/bash

: '
FCC 晶体 {111}<112> 滑移系中的 ISF（本征堆垛层错）, ESF（非本征堆垛层错） 和 TWIN（孪晶层错）模型构建

reference:
- https://mp.weixin.qq.com/s/o0ldM-87tGulOEIBFXSsSw
- https://mp.weixin.qq.com/s/5Mk37UB0SwwbLAw7FaZRZg
'


#-------------------------------- Generate stacking fault configurations --------------------------------
sf_generation() {
  # 设置默认值
  element="${element:-Al}"
  a="${a:-4.041}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  dup_z="${dup_z:-6}"
  vacuum="${vacuum:-15.0}"
  vacuum2=$((2 * vacuum))

  # z 方向的原子层间距
  gap_z=$(echo "scale=15; ${a}*sqrt(3)/3" | bc -l)
  # z 方向长度
  len_z=$(echo "scale=15; ${a}*sqrt(3)*${dup_z}" | bc -l)

  # z 方向长度一半
  len_i=$(echo "scale=15; ${len_z}/2" | bc -l)
  # z 方向长度一半 + 1 * z 方向原子层间距
  len_e=$(echo "scale=15; ${len_i}+${gap_z}" | bc -l)
  # z 方向长度一半 + 2 * z 方向原子层间距
  len_t=$(echo "scale=15; ${len_e}+${gap_z}" | bc -l)

  # a*[112]/6 分 40 次滑移
  ndisp=20
  # 沿 [11-2] 方向每次移动的距离
  disp=$(echo "scale=15; ${a}/sqrt(6)/${ndisp}" | bc -l)


  #-------------------------     创建存储结构文件夹     ---------------------------
  folders=("ISF" "ESF" "TWIN")

  for folder in "${folders[@]}"; do
      if [ -d "${folder}" ]; then
          rm -rf "${folder}"
          mkdir "${folder}"
      else
          mkdir "${folder}"
      fi
  done


  #------------------------- 初始位向超胞构建 ---------------------------
  atomsk --create fcc ${a} ${element} \
         orient "[1-10]" "[11-2]" "[111]" \
         -duplicate ${dup_x} ${dup_y} ${dup_z} \
         -shift 0 0 ${vacuum} \
         -cell add ${vacuum2} z \
         -wrap ${element}_init.lmp


  #------------------------- ISF 模型 ---------------------------
  atomsk ${element}_init.lmp ${element}_0_isf.lmp
  for i in $(seq 1 ${ndisp}); do
    atomsk ${element}_$((i-1))_isf.lmp \
           -shift above ${len_i} z 0 ${disp} 0 \
           -wrap ${element}_${i}_isf.lmp
  done


  #------------------------- ESF 模型 ---------------------------
  atomsk ${element}_${ndisp}_isf.lmp ${element}_0_esf.lmp
  for i in $(seq 1 ${ndisp}); do
    atomsk ${element}_$((i-1))_esf.lmp \
           -shift above ${len_e} z 0 ${disp} 0 \
           -wrap ${element}_${i}_esf.lmp
  done


  #------------------------- TWIN 模型 ---------------------------
  atomsk ${element}_${ndisp}_esf.lmp ${element}_0_twin.lmp
  for i in $(seq 1 ${ndisp}); do
    atomsk ${element}_$((i-1))_twin.lmp \
           -shift above ${len_t} z 0 ${disp} 0 \
           -wrap ${element}_${i}_twin.lmp
  done


  # 移动构型文件至对应的文件夹
  mv ${element}_*_isf.lmp  ISF
  mv ${element}_*_esf.lmp  ESF
  mv ${element}_*_twin.lmp TWIN
}


#-------------------------------- Get help --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-e symbol] [-lc lattice_constant] [-d dup_x dup_y dup_z] [-vac vacuum] [-ni num_interval]"

  echo -e "\nGenerate stacking fault configurations of BCC {110}<111> slip system for VASP calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -e element                 element symbol (default: Ti)"
  echo "    -lc lattice_constant       lattice constant (default: 3.252)"
  echo "    -d dup_x dup_y dup_z       duplicate system in the three directions (default: 1 1 6)"
  echo "    -vac vacuum                vacuum thickness (default: 15.0)"
  echo "    -ni num_interval           number of intervals (default: 10)"

  echo -e "\nExamples:"
  echo "    Default settings: ${script_name}"
  echo "    For pure system: ${script_name} -e Ti -lc 3.252 -d 1 1 6 -vac 15.0 -ni 10"
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
