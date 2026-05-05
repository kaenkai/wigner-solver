set terminal pngcairo enhanced font "Times New Roman,18.0" size 1000,1000
set output "test.png"
set datafile separator ','
set key autotitle columnhead
set multiplot

file = '../test.csv'

##### Constants ######
h = 6.626070040e-34
e = 1.6021766208e-19
c = 299792458

###### Plot parameters ######

set style line 1 lc rgb '#000000' lt 1 lw 2 pt 7 ps 1 dashtype 1 # black
set style line 2 lc rgb '#0000ff' lt 1 lw 2 pt 2 ps 1 dashtype 1 # blue
set style line 3 lc rgb '#FF0000' lt 1 lw 2 pt 4 ps 1 dashtype 1 # red
set style line 4 lc rgb '#006400' lt 1 lw 2 pt 6 ps 1 dashtype 1 # dark-green
set style line 5 lc rgb '#a52a2a' lt 1 lw 2 pt 8 ps 1 dashtype 1 # brown
set style line 6 lc rgb '#9400d3' lt 1 lw 2 pt 3 ps 1 dashtype 1 # purple
set style line 7 lc rgb '#008b8b' lt 1 lw 2 pt 5 ps 1 dashtype 1 # dark cyan
set style line 8 lc rgb '#ffd700' lt 1 lw 2 pt 7 ps 1 dashtype 1 # gold
set style line 9 lc rgb '#FF00ff' lt 1 lw 2 pt 7 ps 1 dashtype 1 # magenta
set style line 10 lc rgb '#000080' lt 1 lw 2 pt 7 ps 1 dashtype 1 # navy

set style line 21 lc rgb '#0000ff' lt 1 lw 2 pt 2 ps 1 dashtype 2 # blue
set style line 31 lc rgb '#FF0000' lt 1 lw 2 pt 4 ps 1 dashtype 2 # red

set style line 22 lc rgb '#0000ff' lt 1 lw 2 pt 2 ps 1 dashtype 4 # blue
set style line 32 lc rgb '#FF0000' lt 1 lw 2 pt 4 ps 1 dashtype 4 # red

set ylabel "y_1"
set y2label "y_2" textcolor "#0000ff"

set mytics
set my2tics
set mxtics

set ytics nomirror
set y2tics nomirror textcolor "#0000ff"

set xrange[:]
set yrange[:]

set format y "% .2E"
set format y2 "% .3f"

#
##### Plot 1 #####
#

set origin 0,0.54
set size 1,0.46

set format x ""
set xlabel ""
# set ylabel "n(x) [cm^{-3}]"
# set y2label "U(x) [eV]" textcolor "#0000ff"
# set logscale y
set key over

plot [:] [:] file u 1:2 w l ls 2 axes x1y2,\
    '' u 1:3 w l ls 1

#
##### Plot 2 #####
#

set origin 0,0
set size 1,0.54

set format y "% .2E"
set format y2 "% .2E"
set format x "%.0f"
set xlabel "x [nm]"
set logscale y
# set logscale y2
set key over

x0 = 35; s = 20; U0 = 0.03
plot [:] file u 1:($6) w l ls 2 axes x1y2,\
    '' u 1:4 w l ls 1,\
    '' u 1:5 w l ls 31
    # '' u 1:7 w l ls 3 axes x1y2,\
    # '' u 1:5 w l ls 1
    # U0*exp(-(x-x0)*(x-x0)/s/s)*(-8./s/s/s/s/s/s*(x-x0)*(x-x0)*(x-x0)+12/s/s/s/s*(x-x0))
    # U0*exp(-(x-x0)*(x-x0)/s/s) w l ls 3
