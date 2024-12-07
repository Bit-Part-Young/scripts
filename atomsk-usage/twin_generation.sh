#!/bin/bash

# 孪晶模型构建
# reference: https://mp.weixin.qq.com/s/018ZIXhlSGUtYBkk5XUN4A

# BCC 孪晶
# TODO: 待验证及确定
# atomsk --create bcc 3.165 W orient "[-211]" "[0-11]" "[111]" -duplicate 5 10 15 W.lmp

# atomsk W.lmp -mirror 0 x -wrap W_mirror.lmp

# atomsk --merge x 2 W_mirror.lmp W.lmp W_twin.lmp

# rm W.lmp W_mirror.lmp


# FCC 孪晶
atomsk --create fcc 4.05 Al orient "[11-2]" "[111]" "[-110]" -duplicate 8 4 4 Al.lmp

atomsk Al.lmp -mirror 0 y -wrap Al_mirror.lmp

atomsk --merge Y 2 Al.lmp Al_mirror.lmp Al_twin.lmp

rm Al.lmp Al_mirror.lmp
