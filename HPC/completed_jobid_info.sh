#!/bin/bash


#-------------------------------- 获取已完成任务 JobId 信息 --------------------------------
completed_jobid_info() {

  day_forward=$1
  # 当前时间
  endtime=$(date +%Y-%m-%d)
  # 当前时间前 N 天
  starttime=$(date -d "-${day_forward} day" +%Y-%m-%d)

  sacct --starttime=${starttime} --endtime=${endtime} --state=COMPLETED --format="JobID%30,JobName,State,Start,End,Elapsed,Workdir%-100" > completed_jobid_info.txt

  # 删除空行
  sed -i '/^$/d' completed_jobid_info.txt
  info=$(grep -v -E 'batch|extern|hydra|\.0|---' completed_jobid_info.txt | sed 's|/dssg/home/acct-mseklt/mseklt-ygj|~|')

  # 前 6 列右对齐，第 7 列左对齐
  echo "${info}" | awk '{printf "%10s %10s %10s %21s %21s %11s  %-50s\n", $1, $2, $3, $4, $5, $6, $7}'

  rm completed_jobid_info.txt

}


#-------------------------------- 获取帮助 --------------------------------
get_help() {
  script_name=$(basename "${0}")

  echo -e "\nUsage: ${script_name} days_forward"

  echo -e "Show completed jobs information in Slurm queue system."

  echo -e "\nOptions:"
  echo "    -h, --help       show this help message and exit"
  echo "    days_forward     days forward"

  echo -e "\nAuthor: SLY."
}


#-------------------------------- 主函数 --------------------------------
if [[ $1 = "-h" || $1 = "--help" ]]; then
  get_help
else
  completed_jobid_info $1
fi
