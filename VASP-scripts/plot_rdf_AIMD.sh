#!/bin/bash

# 获取 AIMD 的 RDF 数据并绘制

# 获取 RDF.dat
echo -e "72\n725\n" | vaspkit > /dev/null


# 绘制 RDF
cat >> plot_tmp.gnu << EOF
set loadpath "~/scripts/cms-scripts/plots"
load "config.gnu"

set output "rdf_AIMD.png"

unset key

set xlabel "r (Å)"
set ylabel "g(r)"
plot "PCF.dat" using 1:2 w lp ps 1 pt 3

unset multiplot
unset output
EOF

gnuplot plot_tmp.gnu

rm plot_tmp.gnu

echo -e "\nRDF plot saved to rdf_AIMD.png."
