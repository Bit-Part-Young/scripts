#!/bin/bash

# 获取 AIMD 的 RDF 数据并绘制

set -eu

# 使用 vaspkit获取 RDF.dat
if [[ -x $(command -v vaspkit) ]]; then
  echo -e "72\n725\n" | vaspkit > /dev/null
else
  echo -e "vaspkit is not installed, please install it first."
fi


if hostname | grep -q sjtu; then
  config_path="~/yangsl/scripts/cms-scripts/plots"
else
  config_path="~/scripts/cms-scripts/plots"
fi

# 绘制 RDF
cat >> .plot.gnu << EOF
set loadpath "${config_path}"
load "config.gnu"

set output "rdf_AIMD.png"

unset key

set xlabel "r (Å)"
set ylabel "g(r)"
plot "PCF.dat" using 1:2 w lp ps 1 pt 3

unset multiplot
unset output
EOF

gnuplot .plot.gnu

rm .plot.gnu

echo -e "\nRDF plot saved to rdf_AIMD.png."
