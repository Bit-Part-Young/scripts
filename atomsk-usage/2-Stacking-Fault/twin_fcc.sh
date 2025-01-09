#!/bin/bash

: '
FCC 孪晶模型构建

reference: https://mp.weixin.qq.com/s/018ZIXhlSGUtYBkk5XUN4A
'


atomsk --create fcc 4.05 Al orient "[11-2]" "[-110]" "[111]" -duplicate 2 1 4 Al.lmp

atomsk Al.lmp -mirror 0 z -wrap Al_mirror.lmp

atomsk --merge z 2 Al.lmp Al_mirror.lmp Al_twin.lmp

rm Al.lmp Al_mirror.lmp
