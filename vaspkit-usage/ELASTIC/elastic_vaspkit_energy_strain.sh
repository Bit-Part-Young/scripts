#!/bin/bash

# 能量-应变法，使用 vaspkit 生成应变构型及输入文件

set -eu


#-------------------------------- 生成 INCAR 文件 --------------------------------
function generate_incar() {
    local encut="${encut:-500}"
    local ediff="${ediff:-1E-06}"
    local ediffg="${ediffg:--1E-02}"
    local kspacing="${kspacing:-0.15}"
    local ismear="${ismear:-1}"

    cat > INCAR << EOF
Global Parameters
SYSTEM   = static
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
ISMEAR   = ${ismear}
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

    natoms=$(sed -n '7p' POSCAR | awk '{ for(i=1; i<=NF; i++) a+=$i; print a}')

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


#-------------------------------- 能量-应变法，使用 vaspkit 生成应变构型及输入文件 --------------------------------
function elastic_vaspkit_energy_strain() {

    local encut="${encut:-500}"
    local ediff="${ediff:-1E-06}"
    local ediffg="${ediffg:--1E-02}"
    local kspacing="${kspacing:-0.15}"
    local ismear="${ismear:-1}"
    local spacegroup="${spacegroup:-}"

    if ! command -v vaspkit &> /dev/null; then
        echo
        echo "Error: vaspkit could not be found, please check if vaspkit is installed!"
        exit 1
    fi

    if [[ -f "CONTCAR" ]] && [[ ! -f "POSCAR" ]]; then
        mv CONTCAR POSCAR
    fi

    # 生成 INCAR 文件，ISIF=2
    generate_incar ${encut} ${ediff} ${ediffg} ${kspacing} ${ismear}
    # 生成空的 KPOINTS 文件（不用，而是用 INCAR 中的 KSPACING 和 KGAMMA 参数）
    touch KPOINTS
    # 生成 POTCAR 文件
    get_psp2.py
    # 生成 job.slurm 文件
    ncpus=$(get_ncpus)
    slurm_generation.py master -nc ${ncpus}

    # 生成 INPUT.in 文件
    cat > INPUT.in << EOF
1                                                              ! 1 for prep-rocessing, 2 for post-processing
3D                                                             ! 2D for slab, 3D for bulk
7                                                              ! number of strain
-0.015 -0.010 -0.005 0.000 0.005 0.010 0.015                   ! magnitude of strain
EOF

    # 生成 vaspkit.in 文件
    cat > vaspkit.in << EOF
2
201
EOF

    # 是否写入具体的空间群至 SYMMETRY.in 文件
    if [[ ! -z "${spacegroup}" ]]; then
        cat > SYMMETRY.in << EOF
# Read Space Group Number from the SYMMETRY.in file if it exists.
    ${spacegroup}          #Space group number of the input structure
EOF
    fi

    # 拷贝 批量检查、提交由 vaspkit 生成的 能量-应变法 弹性常数计算任务的脚本
    fn="/home/yangsl/scripts/cms-scripts/vaspkit-usage/ELASTIC/submit_batch_vasp_jobs.sh"
    cp ${fn} .

    vaspkit < vaspkit.in

    echo
    bash submit_batch_vasp_jobs.sh check
}


#-------------------------------- Get help --------------------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: ${script_name} [-encut INT] [-ediff FLOAT] [-ediffg FLOAT] [-kspacing FLOAT] [-ismear INT] [-sg INT]"

  echo -e "\nGenerate elastic constant calculation input files (energy-strain method; VASP + vaspkit)."

  echo -e "\nOptions:"
  echo "    -h, --help                 show this help message and exit"
  echo "    -encut INT                 ENCUT (default: 500)"
  echo "    -ediff FLOAT               EDIFF (default: 1E-06)"
  echo "    -ediffg FLOAT              EDIFFG (default: -1E-02)"
  echo "    -kspacing FLOAT            KSPACING (default: 0.15)"
  echo "    -ismear INT                ISMEAR (default: 1)"
  echo "    -sg INT                    add specific space group number (default: empty)"

  echo -e "\nExamples:"
  echo "    Default settings: ${script_name}"
  echo "    For higher accuracy: ${script_name} -encut 600 -ediff 1E-08 -kspacing 0.10"
}


#-------------------------------- Main function --------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
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
    -ismear)
      ismear="$2"
      shift 2
      ;;
    -sg)
      spacegroup="$2"
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

elastic_vaspkit_energy_strain "$@"
