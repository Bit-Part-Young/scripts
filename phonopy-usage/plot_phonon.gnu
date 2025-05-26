set terminal pngcairo size 1920, 1080 font 'Arial, 36'
set output "phonon_gnuplot.png"

set ylabel 'Frequency (THz)'

x1 = 0.24746350
x2 = 0.33495510
x3 = 0.59742970
xmax = 0.81173940
ymin = 0
ymax = 10

set xrange [0:xmax]
set yrange [ymin:ymax]

set xtics ("{/Symbol G}" 0, "X" x1, "K" x2, "{/Symbol G}" x3, "L" xmax)
set ytics 2
unset key

# 添加垂直线
set arrow 1 nohead from x1,ymin to x1,ymax lt 2
set arrow 2 nohead from x2,ymin to x2,ymax lt 2
set arrow 3 nohead from x3,ymin to x3,ymax lt 2

plot 'pho.dat' using 1:($2) w l lw 3
