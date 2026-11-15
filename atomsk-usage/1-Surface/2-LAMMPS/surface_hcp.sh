#!/bin/bash

set -eu

# HCP Basal / (0001), Prismatic / (-1100), Pyramidal I / (01-11), Pyramidal II / (11-22) 表面模型构建
# （atomsk，输出 LAMMPS data）
#
# - 取向均为右手系，与 structure-scripts/build_surface.py 的取向表一致（单胞基矢、原子数相同）
# - -d 与 atomsk -dup 一致：nz 为第三基矢周期数，原子层数 = 每周期层数 × nz
# - 真空层：-cell add 沿第三基矢（锥面为倾斜轴）伸长，使法向高度增加 vac（实际间隙约为 vac + 1 个层间距）
# - 终止面为 atomsk 默认（柱面、锥面 I 在宽层间隙处切开）；需选择终止面请用 build_surface.py --term
#
#   表面          orient (x y z)                 每周期层数  每周期原子数  |a| x |b| (Å, Ti a=2.928 c=4.640)
#   basal         [11-20]  [-1100] [0001]         2           4            2.928 x 5.071
#   prismatic     [-1-120] [0001]  [-1100]        4           4            2.928 x 4.640
#   pyramidalI    [-2110]  [1-213] [01-11]        6           6            2.928 x 5.487 (gamma=105.48)
#   pyramidalII   [-1-123] [1-100] [11-22]        10          20           5.487 x 5.071
#
# 常用命令行示例见 ../1-VASP/surface_hcp.sh 头部注释（-d 含义相同）

separator() {
  printf '\n%*s\n\n' "${COLUMNS:-80}" '' | tr ' ' '#'
}

surface_hcp() {
  # 设置默认值 表面类型 元素符号 晶格常数 超胞尺寸 真空层厚度
  surface_type="${surface_type:-basal}"
  symbol="${symbol:-Ti}"
  a="${a:-2.928}"
  c="${c:-4.640}"
  dup_x="${dup_x:-10}"
  dup_y="${dup_y:-10}"
  dup_z="${dup_z:-20}"
  vacuum="${vacuum:-40.0}"
  vacuum_half=$(awk "BEGIN {print ${vacuum}/2}")

  case "${surface_type}" in
    basal) orient=("[11-20]" "[-1100]" "[0001]") ;;          # 单胞 4 个原子，2 个原子层
    prismatic) orient=("[-1-120]" "[0001]" "[-1100]") ;;     # 单胞 4 个原子，4 个原子层
    pyramidalI) orient=("[-2110]" "[1-213]" "[01-11]") ;;    # 单胞 6 个原子，6 个原子层
    pyramidalII) orient=("[-1-123]" "[1-100]" "[11-22]") ;;  # 单胞 20 个原子，10 个原子层
    *)
      echo "Error: unsupported surface type '${surface_type}' (choose from basal, prismatic, pyramidalI, pyramidalII)" >&2
      exit 1
      ;;
  esac

  # 锥面基矢非正交，atomsk 可能询问，echo n 自动回答
  echo n | atomsk --create hcp "${a}" "${c}" "${symbol}" \
    orient "${orient[@]}" \
    -dup "${dup_x}" "${dup_y}" "${dup_z}" \
    -cell add "${vacuum}" z -shift 0 0 "${vacuum_half}" data.lmp

  separator

  structure_folder="hcp-surfaces"
  mkdir -p "${structure_folder}"
  mv data.lmp "${structure_folder}/${surface_type}.lmp"

  echo "HCP ${symbol} ${surface_type} surface model generated: ${structure_folder}/${surface_type}.lmp"
}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-st STR] [-e STR] [-lc FLOAT FLOAT] [-d N N N] [-vac FLOAT]"

  echo -e "\nGenerate HCP Basal / (0001), Prismatic / (-1100), Pyramidal I / (01-11), Pyramidal II / (11-22) surface model for LAMMPS calculation."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -st STR                    surface type: basal, prismatic, pyramidalI, pyramidalII (default: basal)"
  echo "    -e STR                     element symbol (default: Ti)"
  echo "    -lc FLOAT FLOAT            lattice constant a c (default: 2.928 4.640)"
  echo "    -d N N N                   x y z duplication, z in periods (default: 10 10 20)"
  echo "    -vac FLOAT                 vacuum thickness (default: 40.0)"

  echo -e "\nExamples (>= 24 layers):"
  echo "    ${script_name} -st basal -e Ti -lc 2.928 4.640 -d 10 10 12 -vac 40.0"
  echo "    ${script_name} -st prismatic -e Ti -lc 2.928 4.640 -d 10 10 6 -vac 40.0"
  echo "    ${script_name} -st pyramidalI -e Ti -lc 2.928 4.640 -d 10 10 4 -vac 40.0"
  echo "    ${script_name} -st pyramidalII -e Ti -lc 2.928 4.640 -d 10 10 4 -vac 40.0   # 40 layers"
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

surface_hcp
