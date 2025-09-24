#!/bin/bash

set -eu

# 获取 NEP xyz 文件中的能量、力、位力数据

#-------------------------------- 获取能量、应力、力数据 --------------------------------
get_xyz_efs() {
  xyz_fn="${1:-train.xyz}"
  suffix="${2:-dft}"

  awk 'NF == 1' ${xyz_fn} >natoms_${suffix}.dat
  grep -oP 'energy=\K[^ ]+' ${xyz_fn} >energy_${suffix}_tmp.dat
  paste energy_${suffix}_tmp.dat natoms_${suffix}.dat | awk '{printf "%.10f\n", $1/$2}' >energy_${suffix}.dat

  grep -oP 'virial="\K[^"]+' ${xyz_fn} >virial_${suffix}_tmp.dat
  ncol_virial=$(awk 'END { print NF }' virial_${suffix}_tmp.dat)
  if [[ ${ncol_virial} -eq 6 ]]; then
    cp virial_${suffix}_tmp.dat virial_${suffix}.dat
  elif [[ ${ncol_virial} -eq 9 ]]; then
    awk '{print $1, $5, $9, $2, $6, $3}' virial_${suffix}_tmp.dat >virial_${suffix}.dat
  fi

  awk 'NF != 1' ${xyz_fn} >force_${suffix}_tmp.dat
  grep -v "pbc" force_${suffix}_tmp.dat >force_${suffix}_tmp2.dat
  awk '{for(i=NF-2; i<=NF; i++) printf "%.10f ", $i; print ""}' force_${suffix}_tmp2.dat >force_${suffix}.dat

  rm -f energy_${suffix}_tmp.dat force_${suffix}_tmp*.dat virial_${suffix}_tmp.dat

  echo -e "Energy, Force, Virial, Natoms data saved to *_${suffix}.dat.\n"

  echo -e "Energy count: $(wc -l energy_${suffix}.dat | awk '{print $1}')"
  echo -e "Force count: $(wc -l force_${suffix}.dat | awk '{print $1}')"
  echo -e "Virial count: $(wc -l virial_${suffix}.dat | awk '{print $1}')"
  echo -e "Natoms count: $(wc -l natoms_${suffix}.dat | awk '{print $1}')"
}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [xyz_fn] [suffix]"

  echo -e "\nGet the energy, force, virial and natoms data from NEP xyz configuration file."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    xyz_fn                     xyz filename (default: train.xyz)"
  echo "    suffix                     suffix (default: dft)"

  echo -e "\nExamples:"
  echo "    ${script_name}"
  echo "    ${script_name} train.xyz dft"
}

#-------------------------------- 主函数 --------------------------------
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help
  exit 0
elif [[ $# -eq 0 ]]; then
  get_xyz_efs
elif [[ $# -eq 2 ]]; then
  get_xyz_efs "$@"
fi
