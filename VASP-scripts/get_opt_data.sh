#!/bin/bash
# 获取 VASP 弛豫 (opt) 计算目录的离子步、能量、耗时及收敛状态，表格输出

set -eu
export LC_ALL=C

usage() {
  local s; s=$(basename "$0")
  cat <<EOF

Usage: $s [dirs...]

Get ion steps, energy, time cost and convergence state in VASP relaxation calculation dirs with table format.

Options:
  -h, --help      show this help message and exit
  dirs            VASP calculation dirs ("." = all sub-dirs of current dir)

Examples:
    $s \${PWD}       # 当前目录
    $s test1        # 单个目录
    $s test1 test2  # 多个目录
    $s .            # 当前目录下的所有子目录
EOF
}

# 显示名：相对路径原样；绝对路径（去掉 $PWD 前缀后）取最后 3 级
short_name() {
  local d=${1#./}
  if [[ $d != /* ]]; then echo "$d"
  elif [[ $d == "$PWD" ]]; then basename "$d"
  else d=${d#"$PWD"/}; echo "$d" | awk -F/ '{print (NF > 3 ? $(NF-2)"/"$(NF-1)"/"$NF : $0)}'
  fi
}

# 一次遍历 OUTCAR + OSZICAR 提取所需数据（缺失项输出 "-"）
parse_vasp() {
  awk '
    FILENAME == ARGV[1] {
      if (/executed on/ && start == "") start = $(NF-1) " " $NF
      else if (/NIONS =/)               natoms = $NF
      else if (/^ *NSW *=/)             nsw = $3
      else if (/^ *NELM *=/)            { nelm = $3; sub(/;/, "", nelm) }
      else if (/Elapsed time \(sec\)/)  elapsed = $NF
      else if (/reached required accuracy/) conv = 1
      next
    }
    / F= / {
      nion = $1
      for (i = 1; i < NF; i++) if ($i == "E0=") e0 = $(i+1)
      de_ion = ""
      if (match($0, /d E *= *[-+.0-9Ee]+/)) { de_ion = substr($0, RSTART, RLENGTH); sub(/d E *= */, "", de_ion) }
      next
    }
    /^ *[A-Z][A-Z ]*: *[0-9]/ {   # 电子步: DAV: / RMM: / CG : ...
      s = $0; sub(/^ *[A-Z][A-Z ]*: */, "", s); split(s, a, " ")
      nelec = a[1]; eel = a[2]; deel = a[3]
    }
    function v(x) { return x == "" ? "-" : x }
    END {
      OFS = "\t"
      print v(natoms), v(nsw), v(nelm), v(conv), v(elapsed), v(start), v(nion), v(e0), v(de_ion), v(nelec), v(eel), v(deel)
    }' "$1/OUTCAR" "$1/OSZICAR"
}

# 耗时(s)：已结束取 Elapsed time；未结束取 OUTCAR 最后写入时间 - 开始时间
run_seconds() {  # $1: OUTCAR, $2: elapsed, $3: start (YYYY.MM.DD HH:MM:SS)
  local t0 t1
  if [[ $2 != "-" ]]; then printf '%.0f\n' "$2"
  elif t0=$(date -d "${3//./-}" +%s 2>/dev/null) && t1=$(stat -c %Y "$1" 2>/dev/null); then echo $((t1 - t0))
  else echo "-"
  fi
}

# 表格输出（Folder 列宽自适应）
print_table() {
  awk -F'\t' '
    { row[NR] = $0; if (length($1) > w) w = length($1) }
    function hms(t) { return t == "-" ? "-" : sprintf("%02dh %02dm %02ds", t / 3600, t % 3600 / 60, t % 60) }
    END {
      if (w < 6) w = 6
      fmt = "%-" w "s %7s %9s %15s %11s %12s %13s   %s\n"
      line = sprintf("%" (w + 88) "s", ""); gsub(/ /, "-", line)
      print line; printf fmt, "Folder", "natoms", "Ion_Step", "Energy", "dE", "Energy_pa", "TimeCost", "State"; print line
      for (i = 1; i <= NR; i++) { split(row[i], f, "\t"); printf fmt, f[1], f[2], f[3], f[4], f[5], f[6], hms(f[7]), f[8] }
      print line
    }'
}

get_opt_data() {
  local dir natoms nsw nelm conv elapsed start nion e0 de_ion nelec eel deel state
  for dir in "$@"; do
    [[ -s $dir/OUTCAR && -s $dir/OSZICAR ]] || continue
    IFS=$'\t' read -r natoms nsw nelm conv elapsed start nion e0 de_ion nelec eel deel < <(parse_vasp "$dir")

    if [[ $conv == 1 ]]; then state="Converged"
    elif [[ $elapsed != "-" ]]; then state="Unconverged"   # 已结束但未达到收敛精度（NSW 用尽或中断）
    else state="Running"
    fi
    [[ $nion == "-" ]] && nion=0                          # 尚未完成第一个离子步

    awk -v d="$(short_name "$dir")" -v n="$natoms" -v s="$nion" -v e="$e0" -v de="$de_ion" \
        -v t="$(run_seconds "$dir/OUTCAR" "$elapsed" "$start")" -v st="$state" 'BEGIN {
      OFS = "\t"
      print d, n, s, (e == "-" ? "-" : sprintf("%.6f", e)), (de == "-" ? "-" : sprintf("%.3E", de)),
            (e == "-" || n == "-" ? "-" : sprintf("%.6f", e / n)), t, st }'
  done | print_table
}


#-------------------------------- 主函数 --------------------------------
if [[ $# -eq 0 || $1 == "-h" || $1 == "--help" ]]; then
  usage; exit 0
fi
[[ $# -eq 1 && $1 == "." ]] && set -- ./*/

mapfile -t dirs < <(for d in "$@"; do [[ -d $d ]] && printf '%s\n' "${d%/}"; done | sort -V)
get_opt_data "${dirs[@]}"
