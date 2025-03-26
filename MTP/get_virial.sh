#!/bin/bash

# 获取 OUTCAR、train.cfg 文件中的 Stress/Virial 信息

if [[ -f "OUTCAR" ]]; then
  grep -A2 'FORCE on cell' OUTCAR | awk 'NR <= 3'
  grep -A20 'FORCE on cell' OUTCAR | grep 'Total '
  echo
  grep -A2 'FORCE on cell' OUTCAR | awk 'NR <= 3'
  grep -A20 'FORCE on cell' OUTCAR | grep 'in kB'
  echo
  grep -A2 'FORCE on cell' OUTCAR | awk 'NR <= 3'
  grep -A20 'FORCE on cell' OUTCAR | grep 'external pressure'
fi

if [[ -f "train.cfg" ]]; then
  grep -A1 'PlusStress' train.cfg | awk 'NR == 1 || NR % 3 == 2'
fi
