#!/bin/bash

set -eu

# 统计 xyz 文件中的 group 字段类别及其数量
# reference: https://github.com/Kick-H/For_gpumd/blob/master/NEP_related/Count/count_xyz.sh

# 统计 xyz 文件每帧注释行中某个字段 (key="value" 或 key=value) 的取值及帧数
# 纯 awk 实现, 兼容 macOS (BSD) 与 Linux (GNU), 支持带空格的引号值 (如 group="volume scaling")
count_key() {

  key=$1
  xyz_fn=$2

  awk -v key="${key}" '
    skip > 0 { skip--; next }
    header { header = 0; skip = n
             if (match($0, "(^|[[:space:]])" key "=\"[^\"]*\"")) {
               v = substr($0, RSTART, RLENGTH); sub(/^[[:space:]]*[^=]*="/, "", v); sub(/"$/, "", v)
             } else if (match($0, "(^|[[:space:]])" key "=[^[:space:]]+")) {
               v = substr($0, RSTART, RLENGTH); sub(/^[[:space:]]*[^=]*=/, "", v)
             } else { v = "(none)" }
             if (!(v in cnt)) order[++nv] = v
             cnt[v]++; total++; next }
    /^[[:space:]]*$/ { next }
    { n = $1 + 0; header = 1 }
    END { for (i = 1; i <= nv; i++) printf "%s: %d\n", order[i], cnt[order[i]]
          printf "Total: %d.\n", total }
  ' "${xyz_fn}"
}

count_xyz_group() {

  fxyz=$1

  if [ ! -f "${fxyz}" ]; then
    echo "Error: ${fxyz} not found!"
    exit 1
  fi

  echo -e "\n${fxyz} group count:"
  count_key group "${fxyz}"
}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: $script_name [xyz_fn...]"

  echo -e "\nCount group type in xyz file."

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
  count_xyz_group "${xyz_fn}"
done
