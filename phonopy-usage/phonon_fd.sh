#!/bin/bash

# 生成有限位移法计算声子谱的输入文件（VASP + phonopy）

set -eu


#-------------------------------- 生成 INCAR 文件 --------------------------------
function generate_incar() {
    local encut="${encut:-500}"
    local ediff="${ediff:-1E-06}"
    local ediffg="${ediffg:--1E-02}"
    local kspacing="${kspacing:-0.15}"

    cat > INCAR << EOF
Global Parameters
SYSTEM   = Static
ISTART   = 0
ICHARG   = 2
ISPIN    = 1
LCHARG   = .FALSE.
LWAVE    = .FALSE.
PREC     = Accurate
LREAL    = .FALSE.
ENCUT    = ${encut}

KSPACING = ${kspacing}
KGAMMA   = .TRUE.

Electronic Relaxation
ALGO     = Normal
ISMEAR   = 1
SIGMA    = 0.05
EDIFF    = ${ediff}
NELM     = 300
NELMIN   = 6

Ionic Relxation
NSW     = 0
IBRION  = -1
ISIF    = 2
EDIFFG  = ${ediffg}

NPAR     = 4
EOF
}


#-------------------------------- 根据构型文件设置合适的计算所需的核数 --------------------------------
function get_ncpus() {
    local ncpus=20

    natoms=$(sed -n '7p' SPOSCAR | awk '{ for(i=1; i<=NF; i++) a+=$i; print a}')

    if [[ ${natoms} -lt 10 ]]; then
        ncpus=8
    elif [[ ${natoms} -lt 20 ]]; then
        ncpus=20
    elif [[ ${natoms} -lt 40 ]]; then
        ncpus=40
    else
        ncpus=60
    fi
    echo ${ncpus}
}


#-------------------------------- 生成有限位移法计算声子谱的输入文件（VASP + phonopy --------------------------------
function phonon_finite_displacement() {

    local dup_x="${dup_x:-2}"
    local dup_y="${dup_y:-2}"
    local dup_z="${dup_z:-2}"
    local symprec="${symprec:-0.00001}"
    local encut="${encut:-500}"
    local ediff="${ediff:-1E-06}"
    local ediffg="${ediffg:--1E-02}"
    local kspacing="${kspacing:-0.15}"


    if [[ -f "CONTCAR" ]] && [[ ! -f "POSCAR" ]]; then
        mv CONTCAR POSCAR
    fi

    # 根据对称性获取原胞
    phonopy --symmetry --pa auto --tolerance ${symprec}

    # 备份 POSCAR
    if [[ ! -f "POSCAR.orig" ]]; then
        cp POSCAR POSCAR.orig
    fi

    # 将原胞 PPOSCAR 拷贝为 POSCAR
    cp PPOSCAR POSCAR

    # 扩胞，有限位移
    phonopy -d --dim ${dup_x} ${dup_y} ${dup_z} --pa auto

    ncpus=$(get_ncpus)
    for poscar_disp in `ls POSCAR-*`; do
        index=${poscar_disp#*-}

        folder="disp-${index}"
        if [[ ! -d ${folder} ]]; then
            mkdir ${folder}
        fi

        cp ${poscar_disp} ${folder}/POSCAR

        cd ${folder}

        generate_incar ${encut} ${ediff} ${ediffg} ${kspacing}
        get_psp2.py
        slurm_generation.py master -nc ${ncpus}

        cd ..
    done

    echo
    for i in $(ls -d */ | grep ^disp); do echo -n "${i}: "; cd ${i}; check_vaspi; cd ..; done
}


#-------------------------------- Get help --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-d N N N] [-encut INT] [-ediff FLOAT] [-ediffg FLOAT] [-kspacing FLOAT]"

  echo -e "\nGenerate finite displacement phonon calculation input files (VASP + phonopy)."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -d N N N                   duplicate (default: 2 2 2)"
  echo "    -symprec FLOAT             symmetry precision to identify space group (default: 0.00001)"
  echo "    -encut INT                 ENCUT (default: 500)"
  echo "    -ediff FLOAT               EDIFF (default: 1E-06)"
  echo "    -ediffg FLOAT              EDIFFG (default: -1E-02)"
  echo "    -kspacing FLOAT            KSPACING (default: 0.15)"

  echo -e "\nExamples:"
  echo "    Default settings: ${script_name}"
  echo "    For higher accuracy: ${script_name} -d 2 2 2 -encut 600 -ediff 1E-08 -ediffg -1E-02 -kspacing 0.10"
}


#-------------------------------- Main function --------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d)
      dup_x="$2"
      dup_y="$3"
      dup_z="$4"
      shift 4
      ;;
    -symprec)
      symprec="$2"
      shift 2
      ;;
    -encut)
      encut="$2"
      shift 2
      ;;
    -ediff)
      ediff="$2"
      shift 2
      ;;
    -ediffg)
      ediffg="$2"
      shift 2
      ;;
    -kspacing)
      kspacing="$2"
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

phonon_finite_displacement "$@"
