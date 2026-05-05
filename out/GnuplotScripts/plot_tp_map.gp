set terminal pngcairo enhanced font "Times New Roman,12.0" size 800,800
set output 'tpMap.png'

set multiplot

set style line 1 lc rgb '#000000' lw 1 dashtype 1
set style line 2 lc rgb '#000080' lw 1 dashtype 1
set style line 3 lc rgb '#FF0000' lw 2 dashtype 1
set style line 4 lc rgb '#008000' lw 2 dashtype 1
set style line 5 lc rgb '#800000' lw 2 dashtype 1

set pm3d map interpolate 0,0

set border linecolor rgb "white" lw 1

set ylabel 'p [a.u.]' textcolor "black" # offset -0.5,0
set xlabel '{/Symbol t} [s]' # offset 0,-0.5

set xtics textcolor "black" # 0,0.1,0.4
set ytics textcolor "black" -0.04,0.01,0.04

set mxtics
set mytics

set logscale x

set cbrange [:]
set cbtics textcolor "black"
set cblabel 'f(p) [au]' offset 1

set border linecolor rgb "white" lw 1

unset key

set size 0.8,1.1
set origin .04,0

AU_nm = 0.0529
AU_cm2 = 2.8e-17

splot [1e-13:1e-10] [-0.06:0.06] '../tpMap.out' u 1:2:3
# splot [0:2*pi] [0:2*pi] sin(x)*cos(y)


# set palette rgb 7,5,15; # "traditional pm3d\n(black-blue-red-yellow)"
# set palette rgb 3,11,6; # "green-red-violet"
# set palette rgb 23,28,3; # "ocean (green-blue-white)\ntry also other permutations"
# set palette rgb 21,22,23; # "hot (black-red-yellow-white)"
# set palette rgb 30,31,32; # "color printable on gray\n(black-blue-violet-yellow-white)"
# set palette rgb 33,13,10; # "rainbow (blue-green-yellow-red)"
