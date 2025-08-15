#!/bin/bash

# 获取正在运行的任务 JobId 信息
# sacct 需开启（课题组服务器未开启，超算已开启）


running_jobid_info() {

  tmp_fn="running_jobid_info.txt"

  sacct --state=RUNNING --format="JobID%30,JobName,State,Start,End,Elapsed,Workdir%-100" > "${tmp_fn}"

  # 删除空行
  sed -i '/^$/d' running_jobid_info.txt

  info=$(grep -v -E 'batch|extern|hydra|\.0|---' running_jobid_info.txt | sed 's|/dssg/home/acct-mseklt/mseklt-ygj|~|')

  # 前 6 列右对齐，第 7 列左对齐
  echo "${info}" | awk '{printf "%10s %10s %10s %21s %21s %11s  %-50s\n", $1, $2, $3, $4, $5, $6, $7}'

  rm "${tmp_fn}"

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
  running_jobid_info
fi
