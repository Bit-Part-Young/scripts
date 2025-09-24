#!/bin/bash

set -eu

# 将弛豫 OUTCAR 文件多帧转换为 xyz 文件

outcar2xyz_multipleframes(){
  outcar_fn="$1"
  xyz_fn="$2"
  config_type="$3"
  if [ -z "$outcar_fn" ]; then
    outcar_fn="OUTCAR"
  fi

  if [ -z "$xyz_fn" ]; then
    xyz_fn="NEP-dataset.xyz"
  fi

  if [ -z "$config_type" ]; then
    config_type="outcar2xyz"
  fi

  # 清空或创建文件，不写入任何内容
  > "${xyz_fn}"

  # 每个离子步的开始行
  start_line_array=($(sed -n '/aborting loop because EDIFF is reached/=' "${outcar_fn}"))
  # 每个离子步的结束行
  end_line_array=($(sed -n '/[^ML] energy  without entropy/=' "${outcar_fn}"))

  # 元素对应数目
  ion_number_array=($(grep 'ions per type' "${outcar_fn}" | tail -n 1 | awk -F'=' '{print $2}'))
  # 元素符号
  ion_symbol_array=($(grep 'TITEL' "${outcar_fn}" | awk '{print $4}' | awk -F'_' '{print $1}' | awk '!seen[$0]++'))

  # 原子数
  natoms=$(grep 'number of ions' "${outcar_fn}" | awk '{print $12}')

  for ((i=0; i<${#start_line_array[@]}; i++)); do

    start_line=${start_line_array[i]}
    end_line=${end_line_array[i]}

    # 提取离子步信息所在行到临时文件
    temp_fn="outcar.temporary"
    sed -n "${start_line},${end_line}p" "${outcar_fn}" > "${temp_fn}"

    echo "${natoms}" >> "${xyz_fn}"

    # 点阵
    lattice=$(grep -A 7 "VOLUME and BASIS-vectors are now" "$temp_fn" | tail -n 3 | sed 's/-/ -/g' | awk '{print $1,$2,$3}' | xargs)

    # 能量
    energy=$(grep "energy  without entropy" "${temp_fn}" | tail -n 1 | awk '{print $7}')
    # energy=$(grep "free  energy   TOTEN" "${temp_fn}" | tail -1 | awk '{printf "%.10f\n", $5}')

    # virial 数据
    virial=$(grep -A 20 "FORCE on cell =-STRESS" "${temp_fn}" | grep "Total " | tail -n 1 | awk '{print $2,$5,$7,$5,$3,$6,$7,$6,$4}')

    # 添加多个标签实现
    label_array=("element" "group" "description" "tag")
    if [[ $# -eq 0 || $# -eq 3 ]]; then
      echo "config_type=\"${config_type}\" pbc=\"T T T\" energy=${energy} Lattice=\"${lattice}\" virial=\"${virial}\" Properties=species:S:1:pos:R:3:forces:R:3" >> ${xyz_fn}
    elif [[ $# -gt 3 ]]; then
      echo -n "config_type=\"${config_type}\"" >> ${xyz_fn}
      for(( k=4; k<=$#; k++ )); do
        index=$((k-4))
        echo -n " ${label_array[${index}]}=\"${!k}\"" >> "${xyz_fn}"
      done
      echo " pbc=\"T T T\" energy=${energy} Lattice=\"${lattice}\" virial=\"${virial}\" Properties=species:S:1:pos:R:3:forces:R:3" >> ${xyz_fn}
    fi

    # 输出对应个数的元素符号
    for((j=0;j<${#ion_number_array[*]};j++));do
      printf ''${ion_symbol_array[j]}'%.0s\n' $(seq 1 1 ${ion_number_array[j]}) >> "symbol.temporary"
    done

    # 原子位置、力
    grep -A $((natoms + 1)) "TOTAL-FORCE (eV/Angst)" "${temp_fn}" | tail -n "${natoms}" > "positions_forces.temporary"

    # 将元素符号和原子位置、力合并
    paste "symbol.temporary" "positions_forces.temporary" >> "${xyz_fn}"

    # 删除临时文件
    rm -f *.temporary
  done

  echo " ${#start_line_array[@]} frames ${outcar_fn} converted to ${xyz_fn}."
}


#-------------------------------- 获取帮助 --------------------------------
get_help(){
  script_name=$(basename $0)

  echo -e "\nUsage: $script_name [outcar_fn] [xyz_fn] [config_type] [addtional_labels]"

  echo -e "\nConvert opt OUTCAR multiple frames to xyz file for NEP training."

  echo -e "\nOptions:"
  echo "    -h, --help          show this help message and exit"
  echo "    outcar_fn           OUTCAR filename, default OUTCAR"
  echo "    xyz_fn              xyz filename, default NEP-dataset.xyz"
  echo "    config_type         config type tag, default 'outcar2xyz'"
  echo "    addtional_labels    additional labels, default element group description tag"

  echo -e "\nExamples:"
  echo "    $script_name"
  echo "    $script_name OUTCAR NEP-dataset.xyz 'outcar2xyz'"
  echo "    $script_name OUTCAR NEP-dataset.xyz 'outcar2xyz' 'Nb'"
  echo "    $script_name OUTCAR NEP-dataset.xyz 'outcar2xyz' 'Nb' 'relaxation' 'snapshot 1' 'train'"
  echo "    $script_name OUTCAR NEP-dataset.xyz 'outcar2xyz' 'Nb' 'stacking fault' 'configuration 1' 'train'"
  exit 0
}


#-------------------------------- 主函数 --------------------------------
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help

elif [[ $# -eq 0 ]]; then
  outcar2xyz_multipleframes

elif [[ $# -eq 3 ]]; then
  outcar2xyz_multipleframes "$1" "$2" "$3"

elif [[ $# -gt 3 ]]; then
  outcar2xyz_multipleframes "$1" "$2" "$3" "${@:4}"

fi
