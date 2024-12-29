set terminal pngcairo enhanced font "Times New Roman,18.0" size 1000,1000
set output 'map.png'
set datafile separator ','
set key autotitle columnhead

set multiplot

set style line 1 lc rgb '#000000' lw 1 dashtype 1
set style line 2 lc rgb '#000080' lw 1 dashtype 1
set style line 3 lc rgb '#FF0000' lw 2 dashtype 1
set style line 4 lc rgb '#008000' lw 2 dashtype 1
set style line 5 lc rgb '#800000' lw 2 dashtype 1

set pm3d map interpolate 0,0
# set view 70, 90, 1, 1.5

set border linecolor rgb "white" lw 1

set ylabel 'It. nr' textcolor "black" # offset -0.5,0
set xlabel 'x [{/symbol m}m]' # offset 0,-0.5
# set xlabel 'p [au]' # offset 0,-0.5

set xtics textcolor "black" # 0,0.1,0.4
set ytics textcolor "black" # -0.04,0.01,0.04

set mxtics
# set mytics

set yrange [:]
set xrange [:]

set cbrange [:]
set cbtics textcolor "black"
set cblabel 'u^H [eV]' offset 1
# set cblabel 'rho' offset 1

# set logscale cb

set border linecolor rgb "white" lw 1

unset key

set size 1,1.1
set origin 0,-0.04

AU_nm = 0.0529
AU_cm2 = 2.8e-17
AU_eV = 27.211

splot [:] [:] '../poisson_step.out' u ($2*AU_nm*1e-3):($1):($4) with pm3d
# splot [-0.05:0.05] 'poisson_step.out' u 2:($1*0.1):5 with pm3d


# set palette rgb 7,5,15; # "traditional pm3d\n(black-blue-red-yellow)"
# set palette rgb 3,11,6; # "green-red-violet"
# set palette rgb 23,28,3; # "ocean (green-blue-white)\ntry also other permutations"
# set palette rgb 21,22,23; # "hot (black-red-yellow-white)"
# set palette rgb 30,31,32; # "color printable on gray\n(black-blue-violet-yellow-white)"
# set palette rgb 33,13,10; # "rainbow (blue-green-yellow-red)"
