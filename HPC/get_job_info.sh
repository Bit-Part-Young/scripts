#!/bin/bash

# 获取正在运行的任务 JobId 信息


#-------------------------------- 获取 VASP 数据 --------------------------------
get_job_info() {

  # 输出表头
  echo "|-------------------------------------------------------------------------------------"
  echo "|   JobId   |   JobName   |     StartTime  |  RunTime  |     Path          |"
  # echo "|        Path         |   JobId   | JobName | Partition | NumCPUs   |       StartTime       |  RunTime  |   TimeCost    |"
  echo "|-------------------------------------------------------------------------------------"

  job_array=($(squeue | awk 'NR>1 {print $1}'))
  for jobid in "${job_array[@]}"; do
    scontrol show job "${jobid}" > tmp.txt

    cal_path=$(grep 'WorkDir' tmp.txt | awk -F'=' '{print $2}' | awk '{print $1}')
    # 将 cal_path 中的 $HOME 全路径替换为 ~
    cal_path=${cal_path//$HOME/\~}

    jobname=$(grep 'JobName' tmp.txt | awk -F' ' '{print $2}' | awk -F'=' '{print $2}')
    # partition=$(grep 'Partition' tmp.txt | awk -F'=' '{print $2}' | awk '{print $1}')
    # numcpus=$(grep 'NumCPUs' tmp.txt | awk -F'=' '{print $2}' | awk '{print $1}')
    starttime=$(grep 'StartTime' tmp.txt | awk -F' ' '{print $1}' | awk -F'=' '{print $2}')
    runtime=$(grep 'RunTime' tmp.txt | awk -F' ' '{print $1}' | awk -F'=' '{print $2}')

    printf "| %-8s  |     %-4s    |  %-20s  | %-20s | %-35s |\n" "${jobid}" "${jobname}" "${starttime}" "${runtime}" "${cal_path}"
    echo "|-------------------------------------------------------------------------------------"
  done
}

#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "${0}")

  echo -e "\nUsage: ${script_name} num_cpu"

  echo -e "Show running jobs information in HPC."

  echo -e "\nOptions:"
  echo "    -h, --help       show this help message and exit"

  echo -e "\nAuthor: YSL."
}


if [[ $1 = "-h" || $1 = "--help" ]]; then
  get_help
else
  get_job_info
fi
