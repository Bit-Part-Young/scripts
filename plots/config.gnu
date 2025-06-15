# config.gnu gnuplot 绘图配置文件

set terminal pngcairo size 1000,800 enhanced font "Times New Roman,25"

# 设置调色板
# reference: [gnuplot 科技绘图的调色板 - Jerkwin](https://jerkwin.github.io/2018/08/20/%E7%A7%91%E6%8A%80%E7%BB%98%E5%9B%BE%E7%9A%84%E8%B0%83%E8%89%B2%E6%9D%BF/)
setpal="if(pal eq 'cls'){set colorsequence classic};if(pal ne 'def' && pal ne 'cls'){do for[i=1:words(value(pal))]{set style line i lw 4 lc rgb word(value(pal),i)}}"

fav="#1F77B4 #FF7400 #00A13B #D62728 #984EA3 #A65628 #EE0F84 #7F7F7F #BCBD22 #17BECF"
pal='fav'; @setpal

set border lw 5.0

set xtics nomirror
set ytics nomirror

# unset key

# 设置全局线宽（还可设置 点的大小）
set for [i=1:7] style line i lw 5
set style increment user


# 使用示例
# set loadpath "/home/yangsl/scripts/cms-scripts/plots"
# load "config.gnu"

# set output "test.png"

# set xr [0:5]
# plot for [i=1:7] sin(x)+i*0.1 w l
