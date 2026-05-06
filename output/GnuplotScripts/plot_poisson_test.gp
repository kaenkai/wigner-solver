set terminal pngcairo enhanced font "Times New Roman,18.0" size 1000,1000
set output 'poisson_test.png'
set datafile separator ','
set key autotitle columnhead
set multiplot

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

#
# Plot
#

set ylabel 'U [eV]'
set y2label '{/Symbol r}_{new} [C/cm^{3}]'  textcolor '#FF0000'
set xlabel 'x [nm]'
set ytics nomirror
set y2tics nomirror textcolor '#FF0000'
# set yrange [0:]
# set xrange [0:10]
# set logscale y2
set format y2 "%.1E"

set key over

AU_nm = 0.0529; E0 = 1.602E-19; AU_cm2 = 2.8e-17; AU_cm = 5.29e-9

set label "U_{bias} = 0 meV\np=2E-04\nUDS2" at graph 0.05,0.95 textcolor 'black'

# 1: x, 2: U, 3: Delta_U, 4: nC, 5: Current
plot '../poisson_test.csv' u 1:2 w l ls 1, '' u 1:3 w l ls 2, '' u 1:4 w l ls 4,\
    '' u 1:5 w l ls 3 axes x1y2
    # '' u 1:9 w l ls 5 axes x1y2 , [0:4000] 1.33e+04 w l ls 6 axes x1y2
    # '' u 1:5 w l ls 3 axes x1y2, '' u 1:6 w l ls 5 axes x1y2
