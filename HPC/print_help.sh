#!/bin/bash

# 定义 print_help 函数，供其他 Shell 脚本调用

# 检查 echo 是否有 -e 选项
if test "$(echo -e)" = "-e"; then ECHO=echo; else ECHO="echo -e"; fi

#--------------------------- 定义 print_help 函数 -------------------------------
function print_help {
  script_name=$(basename $0)
  $ECHO "${description}."
  $ECHO -n "\nUsage:\n    ${script_name} [-h, --help]"

  # 如果 usage_desc 变量存在，则使用 usage_desc 变量
  if [[ -n $usage_desc ]]; then
    # 方式 1: Usage 中直接写 usage_desc 变量
    $ECHO -n " [${usage_desc}]"
  else
    # 方式 2: 将 options 的 key 写到 Usage 中
    for key in "${!argv_dict[@]}"; do
      $ECHO -n " [${key}]"
      # capitalize the keys in dict
      # $ECHO -n " [${key^^}]"
    done
  fi

  $ECHO "\n\nOptions:"
  printf "%-22s  show this help message and exit\n" "    -h, --help"

  for key in "${!argv_dict[@]}"; do
    printf "%-22s  %s\n" "    ${key}" "${argv_dict[$key]}"
  done

  $ECHO "\nAuthor: YSL."
  exit 0
}
#-------------------------------------------------------------------------------

#----------------------------- 示例代码 -----------------------------------------
: '
# Define the arguments
# 关联数组(Associative Array) 的元素存储是无序的，所以在打印时会出现无序的情况
declare -A argv_dict=(
    ["poscar_file"]="POSCAR file, supported format: *POSCAR*, *.vasp"
    ["type"]="cord type, supported format: frac, cart"
)
description="POSCAR file coordinates format conversion calling pymatgen"
usage_desc="option independent_lattice_constants scale_dimensions species"

# Parse the arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            print_help
            ;;
        *)
            $ECHO "Invalid argument: ${1}. Use '-h/--help' for help.\n"
            print_help
            ;;
    esac
done
'
