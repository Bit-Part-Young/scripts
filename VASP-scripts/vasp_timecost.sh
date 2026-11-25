#!/bin/bash

# 统计 VASP 计算目录的离子步步数、计算耗时、能量信息
# 用法: vasp_timecost.sh [计算目录，默认 .]
# 耗时: 已结束取 OUTCAR 的 Elapsed time；运行中取 已完成离子步 LOOP+ 之和 + 当前离子步已完成电子步 LOOP 之和
# 注: OUTCAR 有写缓冲，最近几个电子步可能尚未写出

set -eu

outcar_fn="${1:-.}/OUTCAR"
[[ -f ${outcar_fn} ]] || { echo "${outcar_fn} not found!" >&2; exit 1; }

# 取 "real time" 之后的数值，避免 cpu time 过大与标签粘连导致列号错位
awk '
  function real(line) { sub(/.*real time/, "", line); return line + 0 }
  /LOOP\+:/                  { nsteps++; loop += real($0); cur = 0; nelec = 0; next }
  /LOOP:/                    { cur += real($0); nelec++ }
  /Elapsed time \(sec\):/    { elapsed = $NF }
  /energy  without entropy=/ { energy = $NF }
  END {
    running = (elapsed == "")
    t = int((running ? loop + cur : elapsed) + 0.5)
    printf "\n%d ionic steps%s, current energy: %s eV, time cost: %dh %dm %ds%s\n",
      nsteps, (running && nelec ? sprintf(" (+%d electronic steps of step %d)", nelec, nsteps + 1) : ""),
      (energy == "" ? "N/A" : sprintf("%.5f", energy)),
      int(t / 3600), int(t % 3600 / 60), t % 60,
      (running ? " (running, sum of LOOP+ and LOOP)" : "")
  }' "${outcar_fn}"
