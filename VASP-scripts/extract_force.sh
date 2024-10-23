#!/bin/bash

# TODO: 提取力数据失败，待检查
# reference: https://github.com/K4ys4r/VASP_Scripts/tree/master/Extract_Forces

Usage() {
  printf "\nExtract forces and atoms positions from VASP OUTCAR File.\n\n"
  printf "Usage:\n"
  printf "    -outcar : Sets the outcar file name.\n"
  printf "    -atom   : Sets the Atom number default value : all atoms.\n"
  printf "    -iter   : Sets the Iterations number default value : all iterations.\n"
  printf "    -out    : Sets the output file.\n"
  printf "    -h      : Print this message. \n"
  printf "    -infos  : Prints infos\n\n"
  printf "Example:\n    Extract_force.sh -atom=4 -iter=all -outcar=OUTCAR\n"
  printf "\nAuthor: Hilal Balout\n"
  printf "E-mail: hilal_balout@hotmail.com\n"
}

for i in "$@"; do
  case $i in
  -atom=*)
    Atom="${i#*=}"
    shift
    ;;
  -iter=*)
    Iter="${i#*=}"
    shift
    ;;
  -outcar=*)
    Inp="${i#*=}"
    shift
    ;;
  -out=*)
    Out="${i#*=}"
    shift
    ;;
  -h)
    Usage
    exit 0
    ;;
  esac
done

if [[ -z "${Atom}" ]]; then
  Atom=all
fi

if [[ -z "${Iter}" ]]; then
  Iter=all
fi

if [[ -z "${Inp}" ]] || [[ ! -f "${Inp}" ]]; then
  printf "\n\t!!!..Error in OUTCAR File Name!\n\n"
  Usage
  exit 0
else
  N_ions=$(grep NIONS ${Inp} | awk '{printf"%d",$NF}')
fi

awk '/TOTAL-FORCE \(eV\/Angst\)/{
    n++
    getline  # skip one line
    if (n=='${Iter}'){
            printf("\nN_iteration : %5d\n",n)
            printf("%5s %10s %10s %10s %12s %12s %12s\n","atom","X","Y","Z","Fx","Fy","Fz")
            for(i=1;i<='${N_ions}';i++){
                getline
                if (i=='${Atom}'){
                    printf("%5d %10.5f %10.5f %10.5f %12.6f %12.6f %12.6f\n",i,$1,$2,$3,$4,$5,$6)
                } else if ("'${Atom}'"=="all") {
                    printf("%5d %10.5f %10.5f %10.5f %12.6f %12.6f %12.6f\n",i,$1,$2,$3,$4,$5,$6)
                }
            }
    } else if ('${Iter}'=="all"){
        printf("\nN_iteration : %5d\n",n)
        printf("%5s %10s %10s %10s %12s %12s %12s\n","atom","X","Y","Z","Fx","Fy","Fz")
        for(i=1;i<='${N_ions}';i++){
            getline
            if (i=='${Atom}'){
                printf("%5d %10.5f %10.5f %10.5f %12.6f %12.6f %12.6f\n",i,$1,$2,$3,$4,$5,$6)
            } else if ('${Atom}'=="all") {
                printf("%5d %10.5f %10.5f %10.5f %12.6f %12.6f %12.6f\n",i,$1,$2,$3,$4,$5,$6)
            }
        }

    }
}' ${Inp} | tee ${Out}
