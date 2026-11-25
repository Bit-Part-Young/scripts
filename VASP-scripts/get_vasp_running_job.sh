#!/bin/bash
# 获取当前用户正在运行的 VASP 作业的计算数据（scf/opt 表格）
#
# 作业目录取自 squeue 的 WORK_DIR：
#   - WORK_DIR 本身含 INCAR → 即为 VASP 计算目录
#   - 否则（如 lane/队列作业在子目录中依次计算）→ 在其下 3 级以内查找
#     最近 ${ACTIVE_MIN} 分钟内有写入且尚未结束的 OUTCAR 所在目录

set -eu

ACTIVE_MIN=${ACTIVE_MIN:-60}
# 只统计该目录下的作业（共享账户如 SiYuan，只看 ~/yangsl）；不存在则不限制
JOB_ROOT=${JOB_ROOT:-$HOME/yangsl}
if [[ -d $JOB_ROOT ]]; then JOB_ROOT=$(readlink -f "$JOB_ROOT"); else JOB_ROOT=""; fi
script_dir=$(dirname "$(readlink -f "$0")")

usage() {
  local s; s=$(basename "$0")
  cat <<EOF

Usage: $s [auto|scf|opt]

Show data of VASP calculations in the current user's RUNNING slurm jobs.

Options:
  -h, --help   show this help message and exit

Arguments:
  auto   按 OUTCAR 中 NSW 自动分为 scf (NSW=0) 与 opt 两张表（默认）
  scf    全部以 get_scf_data.sh 输出
  opt    全部以 get_opt_data.sh 输出

Env:
  ACTIVE_MIN   非 VASP 目录的作业，查找子目录中最近 N 分钟内更新的 OUTCAR（默认 60）
  JOB_ROOT     只统计 WORK_DIR 在该目录下的作业（默认 ~/yangsl，不存在则不限制）

作业名不限：凡 WORK_DIR（或其子目录）中有 VASP 计算文件的运行作业均会被匹配。
EOF
}

cal_type=${1:-auto}
case $cal_type in
  -h|--help) usage; exit 0 ;;
  auto|scf|opt) ;;
  *) echo "Unknown type: $cal_type" >&2; usage >&2; exit 1 ;;
esac

# 调用同目录下的脚本，不存在时用 PATH 中的
run() { if [[ -x $script_dir/$1 ]]; then "$script_dir/$1" "${@:2}"; else "$@"; fi; }

# 收集运行中的 VASP 计算目录（去重）
mapfile -t dirs < <(
  squeue -u "$USER" -h -t RUNNING -o "%Z" | sort -u | while read -r wd; do
    [[ -d $wd ]] || continue
    wd=$(readlink -f "$wd")
    [[ -z $JOB_ROOT || $wd == "$JOB_ROOT" || $wd == "$JOB_ROOT"/* ]] || continue
    if [[ -f $wd/INCAR ]]; then
      echo "$wd"
    else
      find "$wd" -mindepth 2 -maxdepth 4 -name OUTCAR -mmin "-$ACTIVE_MIN" 2>/dev/null |
        while read -r f; do grep -q "Total CPU time used" "$f" || dirname "$f"; done
    fi
  done | sort -uV
)

if [[ ${#dirs[@]} -eq 0 ]]; then
  echo "No running VASP job found for user $USER."
  exit 0
fi

if [[ $cal_type == auto ]]; then
  scf_dirs=() opt_dirs=()
  for d in "${dirs[@]}"; do
    nsw=$(awk '/^ *NSW *=/ {print $3; exit}' "$d/OUTCAR" 2>/dev/null)
    if [[ ${nsw:-0} -eq 0 ]]; then scf_dirs+=("$d"); else opt_dirs+=("$d"); fi
  done
  [[ ${#scf_dirs[@]} -gt 0 ]] && { echo "== SCF (NSW=0)"; run get_scf_data.sh "${scf_dirs[@]}"; }
  [[ ${#opt_dirs[@]} -gt 0 ]] && { echo "== OPT/MD (NSW>0)"; run get_opt_data.sh "${opt_dirs[@]}"; }
  exit 0
fi

run "get_${cal_type}_data.sh" "${dirs[@]}"
