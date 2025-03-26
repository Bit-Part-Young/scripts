#!/bin/bash


# 将 OUTCAR 转换为 xyz 文件（静态计算，单帧），作为 NEP 的训练集构型

# reference: https://github.com/brucefan1983/GPUMD/blob/master/tools/vasp2xyz/outcar2xyz/singleFrame-outcars2nep-exyz.sh

# [ ] 考虑如何添加 config_type 标签


#-------------------------------- 将 OUTCAR 转换为 xyz 文件 --------------------------------
outcar2xyz_singleframe(){
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

  # 原子数
  natoms=$(grep "number of ions" ${outcar_fn} |awk '{print $12}')
  echo $natoms > $xyz_fn

  # 点阵
  lattice=$(grep -A 7 "VOLUME and BASIS-vectors are now" ${outcar_fn} | tail -3 |sed 's/-/ -/g' | awk '{print $1,$2,$3}' |xargs)

  # 能量
  # energy=$(grep "free  energy   TOTEN" ${outcar_fn} | tail -1 | awk '{printf "%.6f\n", $5}')
  energy=$(grep "energy  without entropy" ${outcar_fn} | tail -1 | awk '{printf "%.6f\n", $7}')

  # virial 数据
  virial=$(grep -A 20 "FORCE on cell =-STRESS" ${outcar_fn} | grep "Total " | tail -1 | awk '{print $2,$5,$7,$5,$3,$6,$7,$6,$4}')

  echo "config_type=\"${config_type}\" pbc=\"T T T\" energy=${energy} Lattice=\"${lattice}\" virial=\"${virial}\" Properties=species:S:1:pos:R:3:forces:R:3" >> $xyz_fn

  # 元素对应数目
  ion_number_array=($(grep "ions per type"  ${outcar_fn} | tail -1 | awk -F"=" '{print $2}'))
  # 元素符号
  ion_symbol_array=($(grep "TITEL" ${outcar_fn} | awk '{print $4}' | awk -F"_" '{print $1}' | awk '!seen[$0]++'))

  # 输出对应个数的元素符号
  for((j=0;j<${#ion_number_array[*]};j++)); do
    printf ''${ion_symbol_array[j]}'%.0s\n' `seq 1 1 ${ion_number_array[j]}` >> symbol.temporary
  done

  # 原子位置、力
  grep -A $(($natoms + 1)) "TOTAL-FORCE (eV/Angst)" ${outcar_fn} | tail -n $natoms > positions_forces.temporary

  # 将元素符号和原子位置、力合并
  paste symbol.temporary positions_forces.temporary >> $xyz_fn

  # 删除临时文件
  rm -f *.temporary

  echo "${outcar_fn} converted to ${xyz_fn}."
}


#-------------------------------- 获取帮助 --------------------------------
get_help(){
  script_name=$(basename $0)

  echo -e "\nUsage: $script_name [outcar_fn] [xyz_fn] [config_type]"

  echo -e "\nConvert OUTCAR to xyz file for NEP training."

  echo -e "\nOptions:"
  echo "    -h, --help    show this help message and exit"
  echo "    outcar_fn     OUTCAR filename, default OUTCAR"
  echo "    xyz_fn        xyz filename, default NEP-dataset.xyz"
  echo "    config_type   config type tag, default 'outcar2xyz'"

  echo -e "\nExamples:"
  echo "    $script_name"
  echo "    $script_name OUTCAR NEP-dataset.xyz 'outcar2xyz'"
  exit 0
}


#-------------------------------- 主函数 --------------------------------
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help

elif [[ $# -eq 3 ]]; then
  outcar2xyz_singleframe "$1" "$2" "$3"

elif [[ $# -eq 0 ]]; then
  outcar2xyz_singleframe

fi
