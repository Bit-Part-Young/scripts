#!/bin/bash

# 统计 xyz 文件中的 element 字段类别及其数量
# reference: https://github.com/Kick-H/For_gpumd/blob/master/NEP_related/Count/count_xyz.sh

set -e
set -u

for fxyz in $(ls *.xyz); do

  echo -e "\n$fxyz count info:"

  count_total=0
  elements=$(grep -o 'element=[^ ]*' ${fxyz} | awk -F= '{print $2}' | sort | uniq)
  for element in ${elements[@]}; do
    echo -e "${element}: \c"
    count=$(grep -E element=${element} ${fxyz} | wc -l)
    echo ${count}
    count_total=$((count_total + count))
  done
  echo "Total: ${count_total}"

done
