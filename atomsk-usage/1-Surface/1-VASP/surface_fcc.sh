#!/bin/bash

set -eu

# FCC 结构 (100), (110), (111) 表面模型构建（atomsk，输出 VASP POSCAR）
#
# - 取向均为右手系，与 structure-scripts/build_surface.py 的取向表一致（单胞基矢、原子数相同）
# - -d 与 atomsk -dup 一致：nz 为第三基矢周期数，原子层数 = 每周期层数 × nz
# - 真空层：-cell add 沿 z 增加 vac，顶层与底层（经周期边界）间实际间隙约为 vac + 1 个层间距
#
#   表面  orient (x y z)              每周期层数  每周期原子数  |a| x |b| (Å, Al a=4.041)
#   100   [01-1] [011]   [100]         2           2            2.857 x 2.857
#   110   [-110] [001]   [110]         2           2            2.857 x 4.041
#   111   [1-10] [11-2]  [111]         3           6            2.857 x 4.949
#
# (1) 面内 1x1，原子层数 >= 12（能被 12 整除时恰为 12 层）：
#   surface_fcc.sh -st 100 -e Al -lc 4.041 -d 1 1 6 -vac 15.0     # 12 层, 12 原子
#   surface_fcc.sh -st 110 -e Al -lc 4.041 -d 1 1 6 -vac 15.0     # 12 层, 12 原子
#   surface_fcc.sh -st 111 -e Al -lc 4.041 -d 1 1 4 -vac 15.0     # 12 层, 24 原子
#
# (2) 原子层数 12，原子数 30~70，面内 x、y 尽量等长：
#   surface_fcc.sh -st 100 -e Al -lc 4.041 -d 2 2 6 -vac 15.0     # 48 原子, 5.71 x 5.71
#   surface_fcc.sh -st 110 -e Al -lc 4.041 -d 2 2 6 -vac 15.0     # 48 原子, 5.71 x 8.08
#   surface_fcc.sh -st 111 -e Al -lc 4.041 -d 2 1 4 -vac 15.0     # 48 原子, 5.71 x 4.95

separator() {
  printf '\n%*s\n\n' "${COLUMNS:-80}" '' | tr ' ' '#'
}

surface_fcc() {
  # 设置默认值 表面类型 元素符号 晶格常数 超胞尺寸 真空层厚度
  surface_type="${surface_type:-100}"
  symbol="${symbol:-Al}"
  a="${a:-4.041}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  dup_z="${dup_z:-6}"
  vacuum="${vacuum:-15.0}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")

  case "${surface_type}" in
    100) orient=("[01-1]" "[011]" "[100]") ;;    # 单胞 2 个原子，2 个原子层；或 [010] [001] [100]（4 原子）
    110) orient=("[-110]" "[001]" "[110]") ;;    # 单胞 2 个原子，2 个原子层
    111) orient=("[1-10]" "[11-2]" "[111]") ;;   # 单胞 6 个原子，3 个原子层
    *)
      echo "Error: unsupported surface type '${surface_type}' (choose from 100, 110, 111)" >&2
      exit 1
      ;;
  esac

  atomsk --create fcc "${a}" "${symbol}" \
    orient "${orient[@]}" \
    -dup "${dup_x}" "${dup_y}" "${dup_z}" \
    -cell add "${vacuum}" z -shift 0 0 "${vacuum_half}" \
    -sort species pack -fractional vasp

  separator

  structure_folder="fcc-surfaces"
  mkdir -p "${structure_folder}"
  mv POSCAR "${structure_folder}/POSCAR.${surface_type}"

  echo "FCC ${symbol} (${surface_type}) surface model generated: ${structure_folder}/POSCAR.${surface_type}"
}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-st STR] [-e STR] [-lc FLOAT] [-d N N N] [-vac FLOAT]"

  echo -e "\nGenerate FCC (100), (110), (111) surface model for VASP calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -st STR                    surface type: 100, 110, 111 (default: 100)"
  echo "    -e STR                     element symbol (default: Al)"
  echo "    -lc FLOAT                  lattice constant (default: 4.041)"
  echo "    -d N N N                   x y z duplication, z in periods (default: 1 1 6)"
  echo "    -vac FLOAT                 vacuum thickness (default: 15.0)"

  echo -e "\nExamples (12 layers):"
  echo "    ${script_name} -st 100 -e Al -lc 4.041 -d 1 1 6 -vac 15.0"
  echo "    ${script_name} -st 110 -e Al -lc 4.041 -d 1 1 6 -vac 15.0"
  echo "    ${script_name} -st 111 -e Al -lc 4.041 -d 1 1 4 -vac 15.0"
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

surface_fcc
