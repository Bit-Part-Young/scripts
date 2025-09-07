#!/bin/bash

# 构建 BCC {112}<111> 孪晶模型

element="${element:-Nb}"
a="${a:-3.31}"

atomsk --create bcc ${a} ${element} orient "[111]" "[-110]" "[11-2]" -duplicate 1 1 2 ${element}_cell.xsf

atomsk ${element}_cell.xsf -mirror 0 Z -wrap ${element}_mirror.xsf

atomsk --merge Z 2 ${element}_cell.xsf ${element}_mirror.xsf ${element}_final.cfg

atomsk ${element}_final.cfg vasp

rm *.xsf *.cfg

mv POSCAR ${element}_twin.vasp
