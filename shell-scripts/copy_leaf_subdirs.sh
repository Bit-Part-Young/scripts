#!/bin/bash

# 将文件拷贝到指定目录的所有叶子目录中

set -euo pipefail


copy_leaf_subdirs() {
  local target_dir="$1"
  local source_fn="$2"
  local target_fn="${3:-}"
  local dest_basename
  local count=0
  local dir

  if [[ -n ${target_fn} ]]; then
      dest_basename=${target_fn}
  else
      dest_basename=$(basename ${source_fn})
  fi

  while IFS= read -r -d '' dir; do
      # 跳过包含子目录的目录
      if find "${dir}" -mindepth 1 -maxdepth 1 -type d -print -quit 2>/dev/null | grep -q .; then
          continue
      fi
      cp -- ${source_fn} "${dir}/${dest_basename}"
      count=$((count + 1))
  done < <(find ${target_dir} -type d -print0)

  echo -e "\n${count} files copied to leaf subdirs."
}


#-------------------------------- Help --------------------------------
get_help() {
  script_name=$(basename $0)

  echo -e "\nUsage: ${script_name} <target_dir> <source_fn> [target_fn]"

  echo -e "\nCopy a file into all leaf subdirs."

  echo -e "\nOptions:"
  echo "  -h, --help       show this help message and exit"
  echo "  target_dir       target directory to search"
  echo "  source_fn        source filename"
  echo "  target_fn        target filename (optional)"
}


if [[ $# -ge 1 && ( $1 == "-h" || $1 == "--help" ) ]]; then
  get_help
  exit 0
elif [[ $# -lt 2 || $# -gt 3 ]]; then
  get_help
  exit 1
elif [[ ! -d $1 ]]; then
  echo -e "\nError: not a directory: $1"
  exit 1
elif [[ ! -f $2 ]]; then
  echo -e "\nError: not a regular file: $2"
  exit 1
fi

copy_leaf_subdirs "$1" "$2" "${3:-}"
