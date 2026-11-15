#!/bin/bash

set -eu

# BCC 结构 (100), (110), (111), (112), (123) 表面模型构建（atomsk，输出 VASP POSCAR）
#
# - 取向均为右手系，与 structure-scripts/build_surface.py 的取向表一致（单胞基矢、原子数相同）
# - -d 与 atomsk -dup 一致：nz 为第三基矢周期数，原子层数 = 每周期层数 × nz
# - 真空层：-cell add 沿 z 增加 vac，顶层与底层（经周期边界）间实际间隙约为 vac + 1 个层间距
#
#   表面  orient (x y z)              每周期层数  每周期原子数  |a| x |b| (Å, Nb a=3.307)
#   100   [010]  [001]   [100]         2           2            3.307 x 3.307
#   110   [-110] [001]   [110]         2           4            4.677 x 3.307
#   111   [1-10] [11-2]  [111]         3           6            4.677 x 8.100
#   112   [11-1] [-110]  [112]         6           6            2.864 x 4.677
#   123   [11-1] [-54-1] [123]         14          42           2.864 x 21.432
#
# (1) 面内 1x1，原子层数 >= 12（能被 12 整除时恰为 12 层）：
#   surface_bcc.sh -st 100 -e Nb -lc 3.307 -d 1 1 6 -vac 15.0     # 12 层, 12 原子
#   surface_bcc.sh -st 110 -e Nb -lc 3.307 -d 1 1 6 -vac 15.0     # 12 层, 24 原子
#   surface_bcc.sh -st 111 -e Nb -lc 3.307 -d 1 1 4 -vac 15.0     # 12 层, 24 原子
#   surface_bcc.sh -st 112 -e Nb -lc 3.307 -d 1 1 2 -vac 15.0     # 12 层, 12 原子
#   surface_bcc.sh -st 123 -e Nb -lc 3.307 -d 1 1 1 -vac 15.0     # 14 层, 42 原子
#
# (2) 原子层数 12（123 为 14 层），原子数 30~70，面内 x、y 尽量等长：
#   surface_bcc.sh -st 100 -e Nb -lc 3.307 -d 2 2 6 -vac 15.0     # 48 原子, 6.61 x 6.61
#   surface_bcc.sh -st 110 -e Nb -lc 3.307 -d 1 2 6 -vac 15.0     # 48 原子, 4.68 x 6.61
#   surface_bcc.sh -st 111 -e Nb -lc 3.307 -d 2 1 4 -vac 15.0     # 48 原子, 9.35 x 8.10
#   surface_bcc.sh -st 112 -e Nb -lc 3.307 -d 2 2 2 -vac 15.0     # 48 原子, 5.73 x 9.35
#   surface_bcc.sh -st 123 -e Nb -lc 3.307 -d 1 1 1 -vac 15.0     # 42 原子, 2.86 x 21.43

separator() {
  printf '\n%*s\n\n' "${COLUMNS:-80}" '' | tr ' ' '#'
}

surface_bcc() {
  # 设置默认值 表面类型 元素符号 晶格常数 超胞尺寸 真空层厚度
  surface_type="${surface_type:-100}"
  symbol="${symbol:-Nb}"
  a="${a:-3.307}"
  dup_x="${dup_x:-1}"
  dup_y="${dup_y:-1}"
  dup_z="${dup_z:-6}"
  vacuum="${vacuum:-15.0}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")

  case "${surface_type}" in
    100) orient=("[010]" "[001]" "[100]") ;;     # 单胞 2 个原子，2 个原子层
    110) orient=("[-110]" "[001]" "[110]") ;;    # 单胞 4 个原子，2 个原子层
    111) orient=("[1-10]" "[11-2]" "[111]") ;;   # 单胞 6 个原子，3 个原子层
    112) orient=("[11-1]" "[-110]" "[112]") ;;   # 单胞 6 个原子，6 个原子层
    123) orient=("[11-1]" "[-54-1]" "[123]") ;;  # 单胞 42 个原子，14 个原子层
    *)
      echo "Error: unsupported surface type '${surface_type}' (choose from 100, 110, 111, 112, 123)" >&2
      exit 1
      ;;
  esac

  atomsk --create bcc "${a}" "${symbol}" \
    orient "${orient[@]}" \
    -dup "${dup_x}" "${dup_y}" "${dup_z}" \
    -cell add "${vacuum}" z -shift 0 0 "${vacuum_half}" \
    -sort species pack -fractional vasp

  separator

  structure_folder="bcc-surfaces"
  mkdir -p "${structure_folder}"
  mv POSCAR "${structure_folder}/POSCAR.${surface_type}"

  echo "BCC ${symbol} (${surface_type}) surface model generated: ${structure_folder}/POSCAR.${surface_type}"
}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-st STR] [-e STR] [-lc FLOAT] [-d N N N] [-vac FLOAT]"

  echo -e "\nGenerate BCC (100), (110), (111), (112), (123) surface model for VASP calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -st STR                    surface type: 100, 110, 111, 112, 123 (default: 100)"
  echo "    -e STR                     element symbol (default: Nb)"
  echo "    -lc FLOAT                  lattice constant (default: 3.307)"
  echo "    -d N N N                   x y z duplication, z in periods (default: 1 1 6)"
  echo "    -vac FLOAT                 vacuum thickness (default: 15.0)"

  echo -e "\nExamples (12 layers):"
  echo "    ${script_name} -st 100 -e Nb -lc 3.307 -d 1 1 6 -vac 15.0"
  echo "    ${script_name} -st 110 -e Nb -lc 3.307 -d 1 1 6 -vac 15.0"
  echo "    ${script_name} -st 111 -e Nb -lc 3.307 -d 1 1 4 -vac 15.0"
  echo "    ${script_name} -st 112 -e Nb -lc 3.307 -d 1 1 2 -vac 15.0"
  echo "    ${script_name} -st 123 -e Nb -lc 3.307 -d 1 1 1 -vac 15.0   # 14 layers"
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

surface_bcc
