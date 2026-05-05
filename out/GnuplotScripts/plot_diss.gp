set terminal pngcairo dashed enhanced size 800, 500
set output 'diss.png'

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
# Plot 1
#

set size .99, 1       # set size of inset
set origin 0, 0

set ylabel 'R({/Symbol t})/R({/Symbol t}=0)'  # 'J({/Symbol t})/J_0'
set xlabel '{/Symbol t} [s]'

set key top right

set logscale x

AU_nm = 0.0529; E0 = 1.602E-19; AU_Acm2 = 2.364e14; AU_cm2 = 2.8e-17

file1 = '../diss_tau-tR.out'
file2 = '../diss_tau-tM.out'
file3 = '../diss_tau-tRM.out'

plot file1 u 2:5 w l ls 2 title '{/Symbol t} = {/Symbol t}_R',\
    file2 u 2:5 w l ls 3 title '{/Symbol t} = {/Symbol t}_M',\
    file3 u 2:5 w l ls 4 title '{/Symbol t} = {/Symbol t}_R + {/Symbol t}_M'

#
# Inset
#

set size 0.6, 0.6       # set size of inset
set origin 0.3, 0.2

set ylabel 'J({/Symbol t})/J({/Symbol t}=0)'
unset key
set label "U_{Bias} = 1 meV" at 2e-13,0.9

plot file1 u 2:3 w l ls 2 notitle,\
    file2 u 2:3 w l ls 3 notitle,\
    file3 u 2:3 w l ls 4 notitle
