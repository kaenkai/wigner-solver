set terminal pngcairo enhanced font "Times New Roman,18.0" size 1000,1000
set output 'df_map.png'
# set terminal postscript eps size 6cm,5cm enhanced color # font 'Times New Roman, 12'
# set output 'df_map.eps'

set multiplot

set style line 1 lc rgb '#000000' lw 1 dashtype 1
set style line 2 lc rgb '#000080' lw 1 dashtype 1
set style line 3 lc rgb '#FF0000' lw 2 dashtype 1
set style line 4 lc rgb '#008000' lw 2 dashtype 1
set style line 5 lc rgb '#800000' lw 2 dashtype 1

# set palette rgb 30,31,32
# set hidden3d
# set ticslevel 0.8
# set isosample 50
set pm3d map interpolate 0,0

set border linecolor rgb "white" lw 1
# set grid linecolor "white" lw 2

# set title '{/Symbol t}_R = 1x10^{-12} s, J = 3.71568x10^{06} A/cm^2'

set ylabel 'p [a.u.]' textcolor "black" # offset -0.5,0
set xlabel 'x [{/Symbol m}m]' # offset 0,-0.5

set xtics textcolor "black"
set ytics textcolor "black" # -0.03,0.01,0.03

set mxtics
set mytics

set cbrange [:]
set cbtics textcolor "black"
# set cblabel 'f(x,p;{/Symbol t})-f(x,p;{/Symbol t}=0) [a.u.]' offset 1
set cblabel 'f(x,p) [a.u.]' offset 1

# set format cb "%.1tx10^{%T}"

set border linecolor rgb "white" lw 1
# set grid linecolor "white" lw 2

unset key

set size 1,1.1
set origin 0,-0.04

AU_nm = 0.0529
AU_cm2 = 2.8e-17

# set label "U_{bias} = 100 meV" at graph 0.1,0.9 textcolor "white" front
# set label "{/Symbol t} = {/Symbol t}_M" at graph 0.1,0.85 textcolor "white" front

splot [:] [:] '../wf.out' u ($1*AU_nm/1e3):2:3 with pm3d
# splot [0:2*pi] [0:2*pi] sin(x)*cos(y)

# set palette rgb 7,5,15; # "traditional pm3d\n(black-blue-red-yellow)"
# set palette rgb 3,11,6; # "green-red-violet"
# set palette rgb 23,28,3; # "ocean (green-blue-white)\ntry also other permutations"
# set palette rgb 21,22,23; # "hot (black-red-yellow-white)"
# set palette rgb 30,31,32; # "color printable on gray\n(black-blue-violet-yellow-white)"
# set palette rgb 33,13,10; # "rainbow (blue-green-yellow-red)"
