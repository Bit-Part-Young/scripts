#!/usr/bin/bash

# 获取 MTP cfg 文件中的能量、力、应力数据

#-------------------------------- 获取能量、应力、力数据 --------------------------------
get_mtp_efs() {
  cfg_fn="${1:-train.cfg}"
  suffix="${2:-dft}"

  awk '/Size/ {getline; print $1}' ${cfg_fn} > natoms_${suffix}.dat
  awk '/Energy/ {getline; printf "%.10f\n", $1}' ${cfg_fn} > energy_${suffix}_tmp.dat
  paste energy_${suffix}_tmp.dat natoms_${suffix}.dat | awk '{printf "%.10f\n", $1/$2}' > energy_${suffix}.dat

  awk '/PlusStress/ {getline; for(i=1; i<=NF; i++) printf "%.10f ", $i; print ""}' ${cfg_fn} > stress_${suffix}.dat

  awk '/AtomData/,/Energy/' ${cfg_fn} > forces_${suffix}_tmp.dat
  sed -i -e '/Energy/d' -e '/AtomData/d' forces_${suffix}_tmp.dat
  awk '{for(i=NF-2; i<=NF; i++) printf "%.10f ", $i; print ""}' forces_${suffix}_tmp.dat > forces_${suffix}.dat

  rm -f energy_${suffix}_tmp.dat forces_${suffix}_tmp.dat

  echo -e "Energy, Forces, Stress, Natoms data saved to *_${suffix}.dat.\n"

  echo -e "Energy count: $(wc -l energy_${suffix}.dat | awk '{print $1}')"
  echo -e "Forces count: $(wc -l forces_${suffix}.dat | awk '{print $1}')"
  echo -e "Stress count: $(wc -l stress_${suffix}.dat | awk '{print $1}')"
  echo -e "Natoms count: $(wc -l natoms_${suffix}.dat | awk '{print $1}')"
}

#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [cfg_fn] [suffix]"

  echo -e "\nGet the energy, forces and stress data from MTP cfg format file."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    cfg_fn                     MTP cfg filename (default: train.cfg)"
  echo "    suffix                     suffix (default: dft)"

  echo -e "\nExamples:"
  echo "    ${script_name}"
  echo "    ${script_name} train.cfg dft"
}

#-------------------------------- 主函数 --------------------------------
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  get_help
  exit 0
elif [[ $# -eq 0 ]]; then
  get_mtp_efs
elif [[ $# -eq 2 ]]; then
  get_mtp_efs "$@"
fi
