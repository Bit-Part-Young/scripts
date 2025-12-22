#!/bin/bash


#------------------------- 构建 HCP {10-12} 孪晶模型 --------------------------
function twin_10_12_generation() {

  symbol=${symbol:-Mg}
  a="${a:-3.200}"
  c="${c:-5.160}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  dup_z="${dup_z:-2}"
  dup_z2=$((dup_z * 2))

  atomsk --create hcp ${a} ${c} ${symbol} orient "[1-210]" "[-1011]" "[10-12]" -duplicate ${dup_x} ${dup_y} ${dup_z2} -sort z up -wrap -ow vasp

  mv POSCAR bulk_${symbol}.vasp

  atomsk --create hcp ${a} ${c} ${symbol} orient "[1-210]" "[-1011]" "[10-12]" -duplicate ${dup_x} ${dup_y} ${dup_z} -ow ${symbol}_cell.xsf
  atomsk --create hcp ${a} ${c} ${symbol} orient "[1-210]" "[-1011]" "[-101-2]" -duplicate ${dup_x} ${dup_y} ${dup_z} -ow ${symbol}_mirror.xsf
  atomsk --merge Z 2 ${symbol}_mirror.xsf ${symbol}_cell.xsf -ow ${symbol}_final.cfg
  atomsk ${symbol}_final.cfg -sort z up -wrap -ow vasp

  rm *.xsf *.cfg

  mv POSCAR twin_10_12_${symbol}_original.vasp
  if [[ ! -f "twin_10_12_${symbol}_final.vasp" ]]; then
    cp twin_10_12_${symbol}_original.vasp twin_10_12_${symbol}_final.vasp
  fi


  echo -e "\nHCP ${symbol} {10-12} twin model has been generated, saved as bulk_${symbol}.vasp and twin_10_12_${symbol}_original.vasp."
  echo -e "\nNote: You need manually modify the atomic positions in twin_10_12_${symbol}_final.vasp to get the final twin model."


  : '
      # 原子位置修改
      # 注: 将 z 坐标变成 0；y 坐标的原子间距与孪晶下半部分的保持一致（与中间孪晶的 y 坐标的原子间距是不一样的；不用特别精确）
      0.00000000        4.73226959       22.03063675  ->  0.00000000        4.46246000       0.0

      0.00000000        4.73226959       23.28953039  ->  0.00000000        4.73226959       22.66008357

  '
}


#-------------------------------- Get help --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-e STR] [-lc FLOAT FLOAT] [-d N N N]"

  echo -e "\nGenerate HCP {10-12} twin model."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -e STR                     element symbol (default: Mg)"
  echo "    -lc FLOAT FLOAT            lattice constant (default: 3.200 5.160)"
  echo "    -d N N N                   duplicate structure in the three directions (default: 1 1 2)"

  echo -e "\nExamples:"
  echo "    Default settings: ${script_name}"
  echo "    For HCP Mg system: ${script_name} -e Mg -lc 3.200 5.160 -d 1 1 2"
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
      c="$3"
      shift 3
      ;;
    -d)
      dup_x="$2"
      dup_y="$3"
      dup_z="$4"
      shift 4
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
twin_10_12_generation "$@"
