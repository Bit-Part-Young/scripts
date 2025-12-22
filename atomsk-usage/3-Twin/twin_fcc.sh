#!/bin/bash


#------------------------- 构建 FCC {111}<112> 孪晶模型 --------------------------
function twin_fcc() {

  symbol="${symbol:-Al}"
  a="${a:-4.041}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  dup_z="${dup_z:-2}"
  dup_z2=$((dup_z * 2))

  atomsk --create fcc ${a} ${symbol} orient "[11-2]" "[-110]" "[111]" -duplicate ${dup_x} ${dup_y} ${dup_z2} -sort z up -wrap -ow vasp

  mv POSCAR bulk_${symbol}.vasp

  atomsk --create fcc ${a} ${symbol} orient "[11-2]" "[-110]" "[111]" -duplicate ${dup_x} ${dup_y} ${dup_z} -ow ${symbol}_cell.xsf
  atomsk --create fcc ${a} ${symbol} orient "[11-2]" "[-110]" "[-1-1-1]" -duplicate ${dup_x} ${dup_y} ${dup_z} -ow ${symbol}_mirror.xsf

  atomsk --merge Z 2 ${symbol}_cell.xsf ${symbol}_mirror.xsf ${symbol}_final.cfg

  atomsk ${symbol}_final.cfg -sort z up -wrap -ow vasp

  rm *.xsf *.cfg

  mv POSCAR twin_${symbol}.vasp

  echo -e "\nFCC ${symbol} {111}<112> twin model has been generated, saved as bulk_${symbol}.vasp and twin_${symbol}.vasp."
}


#-------------------------------- Get help --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-e STR] [-lc FLOAT FLOAT] [-d N N N]"

  echo -e "\nGenerate FCC {111}<112> twin model."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -e STR                     element symbol (default: Al)"
  echo "    -lc FLOAT FLOAT            lattice constant (default: 4.041)"
  echo "    -d N N N                   duplicate structure in the three directions (default: 1 1 2)"

  echo -e "\nExamples:"
  echo "    Default settings: ${script_name}"
  echo "    For FCC Al system: ${script_name} -e Al -lc 4.041 -d 1 1 2"
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


twin_fcc "$@"
