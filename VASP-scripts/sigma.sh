#!/bin/bash

: '
检查 entropy T*S 是否小于 1 meV/atom
用法: sigma.sh [计算目录，默认 .]

reference: https://github.com/tamaswells/VASP_script/blob/master/sigma.sh
'

set -eu

dir=${1:-.}
for f in OUTCAR INCAR; do
  [[ -f ${dir}/${f} ]] || { echo "${dir}/${f} not found!" >&2; exit 1; }
done

# 只取未被注释的 SIGMA 设置，缺省为 VASP 默认值 0.2
sigma=$(awk -F'=' 'toupper($1) ~ /^[ \t]*SIGMA[ \t]*$/ {split($2, a, /[ \t;!#]+/); v = (a[1] == "" ? a[2] : a[1])} END {print (v == "" ? 0.2 : v)}' "${dir}/INCAR")

awk -v sigma="${sigma}" '
  /entropy T\*S/ { entropy = $NF }
  /NIONS/        { natoms = $NF }
  END {
    if (entropy == "" || natoms == "") { print "entropy T*S or NIONS not found in OUTCAR!"; exit 1 }
    value = entropy / natoms; if (value < 0) value = -value
    if (value < 0.001)
      printf "Your SIGMA %s is okay; for average entropy T*S %.6f is < 0.001 eV/atom.\n", sigma, value
    else
      printf "Your SIGMA %s is BAD; for average entropy T*S %.6f is > 0.001 eV/atom. Try to decrease SIGMA!\n", sigma, value
  }' "${dir}/OUTCAR"
