#!/bin/bash

: '
提取 VASP OUTCAR 文件中指定离子步中的指定原子的位置、受力信息

reference: https://github.com/K4ys4r/VASP_Scripts/tree/master/Extract_Forces

使用:
  # 提取所有离子步中的所有原子的位置、受力信息，并保存到文件中
  extract_force.sh -outcar=OUTCAR -iter=1 -atom=1 -out=output.txt
  extract_force.sh -outcar=OUTCAR -out=output.txt
  # 提取第 1 个离子步中的所有原子的位置、受力信息
  extract_force.sh -outcar=OUTCAR -iter=1
  # 提取第 1 个离子步中的第 1 个 原子的位置、受力信息
  extract_force.sh -outcar=OUTCAR -iter=1 -atom=1
'


#------------------     帮助信息     ------------------
get_help() {
  script_name=$(basename "$0")

  echo -e "\nUsage: $script_name [-outcar=OUTCAR] [-atom=ATOM] [-iter=ITER] [-out=OUTPUT]"

  echo -e "\nExtract forces and atoms positions from VASP OUTCAR File."

  echo -e "\nOptions:"
  echo "  -h, --help      Show this help message and exit"
  echo "  -outcar OUTCAR  Sets the outcar file name"
  echo "  -atom=ATOM      Sets the Atom number default value : all atoms"
  echo "  -iter=ITER      Sets the Iterations number default value : all iterations"
  echo "  -out=OUTPUT     Sets the output file"
  echo "  -infos          Prints infos"

  echo -e "\nExamples:"
  echo "    ${script_name} -atom=4 -iter=all -outcar=OUTCAR"
  echo "    ${script_name} -atom=4 -iter=1 -outcar=OUTCAR -out=output.log"

  echo -e "\nAuthor: Hilal Balout"
  echo "E-mail: hilal_balout@hotmail.com"
}


#------------------     命令行参数解析     ------------------
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
    get_help
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
  echo -e "\n\t!!!..Error in OUTCAR File Name!\n\n"
  get_help
  exit 0
else
  N_ions=$(grep NIONS "${Inp}" | awk '{printf"%d",$NF}')
fi


#------------------     提取原子位置、受力信息     ------------------
# 注意单双引号的使用
awk '/TOTAL-FORCE \(eV\/Angst\)/{
    n++
    getline  # skip one line
    if (n=="'${Iter}'"){
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
    } else if ("'${Iter}'"=="all"){
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

    }
}' "${Inp}" | tee "${Out}"
