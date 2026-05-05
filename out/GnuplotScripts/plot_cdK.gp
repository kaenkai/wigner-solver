set terminal pngcairo dashed enhanced size 800, 500
set output 'cdK.png'

set multiplot

set style line 1 lc rgb '#000000' lw 2 dashtype 1
set style line 2 lc rgb '#000080' lw 2 dashtype 1
set style line 3 lc rgb '#FF0000' lw 2 dashtype 1
set style line 4 lc rgb '#008000' lw 2 dashtype 1
set style line 5 lc rgb '#800000' lw 2 dashtype 1

#
# Plot
#

set ylabel 'f [a.u.]'
set xlabel 'p [a.u.]'

# set yrange [0:]
# set xrange [-0.1:0.1]
# set logscale y

set key top left

AU_nm = 0.0529; E0 = 1.602E-19; AU_cm2 = 2.8e-17; AU_cm = 5.29e-9

plot '../cdK.out' u 1:2 w l ls 1 notitle
