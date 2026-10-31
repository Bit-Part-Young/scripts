#!/bin/bash

# 生成 VASP + phonopy 声子计算的输入文件：有限位移法（fd，每个位移一个 disp-XXX/）或 DFPT（dfpt，一个 DFPT/）
# 需要 phonopy ≥ 4（phonopy-init），以及 cms-scripts 中的 get_psp2.py、slurm_generation.py、check_vaspi

set -eu

mode="fd"
dim=(2 2 2)
encut=500
ediff="1E-06"
ediffg="-1E-02"
kspacing=0.15
ismear=""
symprec=0.00001
platform="master"


#-------------------------------- 生成 INCAR 文件 --------------------------------
function generate_incar() {
  if [[ ${mode} == "fd" ]]; then
    local system="Static" nsw=0 ibrion=-1 extra="NPAR     = 4"
  else
    # DFPT 不支持 NPAR
    local system="Phonon" nsw=1 ibrion=8 extra="ADDGRID  = .TRUE."
  fi

  cat > INCAR << EOF
Global Parameters
SYSTEM   = ${system}
ISTART   = 0
ICHARG   = 2
ISPIN    = 1
LCHARG   = .FALSE.
LWAVE    = .FALSE.
PREC     = Accurate
LREAL    = .FALSE.
ENCUT    = ${encut}
${extra}

KSPACING = ${kspacing}
KGAMMA   = .TRUE.

Electronic Relaxation
ALGO     = Normal
ISMEAR   = ${ismear}
SIGMA    = 0.05
EDIFF    = ${ediff}
NELM     = 300
NELMIN   = 6

Ionic Relaxation
NSW      = ${nsw}
IBRION   = ${ibrion}
ISIF     = 2
EDIFFG   = ${ediffg}
EOF
}


#-------------------------------- 根据超胞原子数设置核数 --------------------------------
function get_ncpus() {
  local natoms
  natoms=$(sed -n '7p' SPOSCAR | awk '{ for(i=1; i<=NF; i++) a+=$i; print a }')

  if [[ ${natoms} -lt 10 ]]; then
    echo 8
  elif [[ ${natoms} -lt 20 ]]; then
    echo 20
  elif [[ ${natoms} -lt 40 ]]; then
    echo 40
  else
    echo 60
  fi
}


#-------------------------------- 准备一个计算目录 --------------------------------
function setup_job() {
  local folder=$1 poscar=$2 ncpus=$3

  mkdir -p "${folder}"
  cp "${poscar}" "${folder}/POSCAR"
  (
    cd "${folder}"
    generate_incar
    get_psp2.py
    slurm_generation.py "${platform}" -nc "${ncpus}"
    echo -n "${folder}: "
    check_vaspi
  )
}


#-------------------------------- 生成声子计算输入文件 --------------------------------
function phonon_vasp() {
  if ! command -v phonopy-init &> /dev/null; then
    echo -e "\nError: phonopy-init could not be found, please install phonopy >= 4!"
    exit 1
  fi

  if [[ -f "CONTCAR" ]] && [[ ! -f "POSCAR" ]]; then
    mv CONTCAR POSCAR
  fi

  # 根据对称性获取原胞，备份 POSCAR 后以原胞 PPOSCAR 代替
  phonopy-init --symmetry --pa auto --tolerance "${symprec}"
  if [[ ! -f "POSCAR.orig" ]]; then
    cp POSCAR POSCAR.orig
  fi
  cp PPOSCAR POSCAR

  # 扩胞（DFPT 只用 SPOSCAR，有限位移另生成 POSCAR-XXX）
  phonopy-init -d --dim "${dim[@]}" --pa auto

  local ncpus
  ncpus=$(get_ncpus)
  echo
  if [[ ${mode} == "fd" ]]; then
    for poscar_disp in POSCAR-*; do
      setup_job "disp-${poscar_disp#POSCAR-}" "${poscar_disp}" "${ncpus}"
    done
  else
    setup_job "DFPT" "SPOSCAR" "${ncpus}"
  fi
}


#-------------------------------- Get help --------------------------------
function get_help() {
  local script_name
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-m fd|dfpt] [-d N N N] [-encut INT] [-ediff FLOAT] [-ediffg FLOAT] [-kspacing FLOAT] [-ismear INT] [-symprec FLOAT] [-platform STR]"

  echo -e "\nGenerate VASP + phonopy phonon calculation input files (finite displacement or DFPT)."

  echo -e "\nOptions:"
  echo "  -h, --help                 show this help message and exit"
  echo "  -m fd|dfpt                 method: finite displacement or DFPT (default: fd)"
  echo "  -d N N N                   supercell size (default: 2 2 2)"
  echo "  -encut INT                 ENCUT (default: 500)"
  echo "  -ediff FLOAT               EDIFF (default: 1E-06)"
  echo "  -ediffg FLOAT              EDIFFG (default: -1E-02)"
  echo "  -kspacing FLOAT            KSPACING (default: 0.15)"
  echo "  -ismear INT                ISMEAR (default: 1 for fd, 0 for dfpt)"
  echo "  -symprec FLOAT             symmetry precision to identify space group (default: 0.00001)"
  echo "  -platform STR              platform (default: master; master or sy)"

  echo -e "\nExamples:"
  echo "    Default settings: ${script_name}"
  echo "    DFPT: ${script_name} -m dfpt -d 3 3 3"
  echo "    For higher accuracy: ${script_name} -d 4 4 4 -encut 600 -ediff 1E-08 -kspacing 0.10"
}


#-------------------------------- Main function --------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -m)        mode="$2";              shift 2 ;;
    -d)        dim=("$2" "$3" "$4");   shift 4 ;;
    -encut)    encut="$2";             shift 2 ;;
    -ediff)    ediff="$2";             shift 2 ;;
    -ediffg)   ediffg="$2";            shift 2 ;;
    -kspacing) kspacing="$2";          shift 2 ;;
    -ismear)   ismear="$2";            shift 2 ;;
    -symprec)  symprec="$2";           shift 2 ;;
    -platform) platform="$2";          shift 2 ;;
    -h | --help) get_help; exit 0 ;;
    *)         get_help; exit 1 ;;
  esac
done

case "${mode}" in
  fd)   ismear="${ismear:-1}" ;;
  dfpt) ismear="${ismear:-0}" ;;
  *)    echo "Error: -m must be fd or dfpt"; exit 1 ;;
esac

phonon_vasp
