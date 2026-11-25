#!/bin/bash

set -eu

# 统计 xyz 文件构型帧数及总原子数
# reference: https://github.com/Kick-H/For_gpumd/blob/master/NEP_related/Count/count_xyz.sh


count_xyz() {

  xyz_fn=$1

  if [ ! -f "${xyz_fn}" ]; then
    echo "Error: ${xyz_fn} not found!"
    exit 1
  fi

  # 逐帧解析: 第一行为原子数, 第二行为注释行, 之后 natoms 行为原子坐标
  # 纯 awk 实现, 兼容 macOS (BSD) 与 Linux (GNU), 不依赖 wc 的输出格式
  awk -v fn="${xyz_fn}" '
    skip > 0 { skip--; next }
    /^[[:space:]]*$/ { next }
    { n = $1 + 0; frames++; atoms += n; skip = n + 1 }
    END { printf "%s: %d frames, %d atoms.\n", fn, frames, atoms }
  ' "${xyz_fn}"

}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [xyz_fn...]"

  echo -e "\nCount frames and atoms in xyz files."

  echo -e "\nOptions:"
  echo "  -h, --help      show this help message and exit"
  echo "  xyz_fn          xyz filename"

  echo -e "\nExamples:"
  echo "    ${script_name}                # 当前目录下的所有 xyz 文件"
  echo "    ${script_name} 1.xyz          # 单个 xyz 文件"
  echo "    ${script_name} 1.xyz 2.xyz    # 多个 xyz 文件"

  exit 0
}


#-------------------------------- 主函数 --------------------------------
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  get_help
fi

# 未指定文件时, 处理当前目录下的所有 xyz 文件
if [ $# -eq 0 ]; then
  shopt -s nullglob
  xyz_files=(*.xyz)
  shopt -u nullglob
  if [ ${#xyz_files[@]} -eq 0 ]; then
    echo "Error: no xyz file found in current directory!"
    exit 1
  fi
  set -- "${xyz_files[@]}"
fi

for xyz_fn in "$@"; do
  count_xyz "${xyz_fn}"
done
