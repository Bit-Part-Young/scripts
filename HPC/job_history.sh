#!/bin/bash

# 查看提交至 Slurm 队列系统的 job 历史
# sacct 需 enabled

start=$1     # eg. "2025-06-29"
end=$2       # eg. "2025-07-02"

sacct --starttime=${start} --endtime=${end} \
      --format=JobID%10,JobName%9,ReqCPUs,State%7,Start,End,Elapsed,Workdir%-150 \
      | grep -viE 'batch|extern|hydra|-----' \
      | sed "s:/dssg.*-ygj:~:g" \
      | sed 's/[[:space:]]\+$//' > temp.job_info

flag0=$(cat temp.job_info | awk '$4 ~ /COMP/' | wc -l)
echo -e "\nCOMPLETED jobs info:\n"
sed -n 1p temp.job_info
cat temp.job_info | awk '$4 ~ /COMP/'

echo
printf "%`tput cols`s" | tr ' ' '#'
echo

flag1=$(cat temp.job_info | awk '$4 ~ /RUNNING/' | wc -l)
echo -e "\nRUNNING jobs info:\n"
sed -n 1p temp.job_info
cat temp.job_info | awk '$4 ~ /RUNNING/' | sed 's/          Unknown/Unknown/g'

echo
printf "%`tput cols`s" | tr ' ' '#'
echo

flag2=$(cat temp.job_info | awk '$4 ~ /FAILED/' | wc -l)
echo -e "\nFAILED jobs info:\n"
sed -n 1p temp.job_info
cat temp.job_info | awk '$4 ~ /FAILED/'

echo
printf "%`tput cols`s" | tr ' ' '#'
echo

echo -e "\nCOMPLETED ${flag0}; RUNNING ${flag1}; FAILED ${flag2}\n"

rm temp.job_info
