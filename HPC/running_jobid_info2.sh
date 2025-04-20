#!/bin/bash

# 针对 sacct 无法使用的情况，使用 squeue

#-------------------------------- 获取正在运行的任务 JobId 信息 --------------------------------
running_jobid_info2() {

  # 输出表头
  printf "%10s %10s %10s %21s %11s  %-50s\n" "JobId" "JobName" "State" "StartTime" "Elapsed" "WorkDir"

  job_array=($(squeue -u yangsl | awk 'NR>1 {print $1}'))
  for jobid in "${job_array[@]}"; do
    tmp_fn="running_jobid_info2.txt"
    scontrol show job "${jobid}" >"${tmp_fn}"

    workdir=$(grep 'WorkDir' "${tmp_fn}" | awk -F'=' '{print $2}' | awk '{print $1}')
    # 将 workdir 中的 $HOME 全路径替换为 ~
    workdir=${workdir//$HOME/\~}

    jobname=$(grep 'JobName' "${tmp_fn}" | awk -F' ' '{print $2}' | awk -F'=' '{print $2}')
    state=$(grep 'JobState' "${tmp_fn}" | awk -F' ' '{print $1}' | awk -F'=' '{print $2}')
    # partition=$(grep 'Partition' tmp.txt | awk -F'=' '{print $2}' | awk '{print $1}')
    # numcpus=$(grep 'NumCPUs' tmp.txt | awk -F'=' '{print $2}' | awk '{print $1}')
    starttime=$(grep 'StartTime' "${tmp_fn}" | awk -F' ' '{print $1}' | awk -F'=' '{print $2}')
    runtime=$(grep 'RunTime' "${tmp_fn}" | awk -F' ' '{print $1}' | awk -F'=' '{print $2}')

    if [[ ${state} = "RUNNING" ]]; then
      printf "%10s %10s %10s %21s %11s  %-50s\n" "${jobid}" "${jobname}" "${state}" "${starttime}" "${runtime}" "${workdir}"
    fi

    rm "${tmp_fn}"
  done
}

#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "${0}")

  echo -e "\nUsage: ${script_name} num_cpu"

  echo -e "Show running jobs information in Slurm queue system."

  echo -e "\nOptions:"
  echo "    -h, --help       show this help message and exit"

  echo -e "\nAuthor: SLY."
}

#-------------------------------- 主函数 --------------------------------
if [[ $1 = "-h" || $1 = "--help" ]]; then
  get_help
else
  running_jobid_info2
fi
