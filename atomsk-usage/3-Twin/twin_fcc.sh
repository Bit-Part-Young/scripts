#!/bin/bash

# 构建 FCC {111}<112> 孪晶模型

element="${element:-Al}"
a="${a:-4.02}"

atomsk --create fcc ${a} ${element} orient "[11-2]" "[-110]" "[111]" -duplicate 1 1 2 ${element}_cell.xsf

atomsk ${element}_cell.xsf -mirror 0 Z -wrap ${element}_mirror.xsf

atomsk --merge Z 2 ${element}_cell.xsf ${element}_mirror.xsf ${element}_final.cfg

atomsk ${element}_final.cfg POSCAR

rm *.xsf *.cfg

mv POSCAR ${element}_twin.vasp
