set terminal pngcairo enhanced font "Times New Roman,18.0" size 1200,1000
set output "trChar.png"

set datafile separator ','
set key autotitle columnhead

set multiplot

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

set xlabel "It. nr"
# set ylabel "Gęstość prądu [Acm^{-2}]"  # Gęstość prądu [Acm^{-2}]
# set ylabel "dU/U"

set mytics
set my2tics
set mxtics

set format y "% .1tx10^{%.1T}"

# set ytics nomirror textcolor '#0000ff   '
# set y2tics nomirror textcolor '#FF0000'
# set xtics 0,8000,40000

set xrange[0:]
set yrange[:]

set key bottom right # title 'Schemat r.'
set grid

# set label "{/symbol a}=1E-07" at graph 0.05,0.9 textcolor 'black'

#
##### Plot 1 #####
#

set origin 0,0
set size 0.5,0.5
i = 2
plot '../tr_char.csv' u ($1*1):i w l ls 2

#
##### Plot 2 #####
#

set origin 0.48,0
set size 0.5,0.5
i = 3
plot '../tr_char.csv' u ($1*1):i w l ls 2

set logscale y

#
##### Plot 3 #####
#

set origin 0.0,0.5
set size 0.5,0.5
i = 7
plot '../tr_char.csv' u ($1*1):i w l ls 2

#
##### Plot 4 #####
#

set origin 0.48,0.5
set size 0.5,0.5
i = 6
plot '../tr_char.csv' u ($1*1):i w l ls 2

# plot '../tr_char_p001.csv' u 1:i w l ls 2 title 'UDS3',\
#     '../tr_char_p001_HDS22.csv' u 1:i w l ls 3 title 'HDS22'
    # '../tr_char_p0008.csv' u 1:i w l ls 2 title '0.008',\
    # '../tr_char_p0006.csv' u 1:i w l ls 3 title '0.006',\
    # '../tr_char_p0004.csv' u 1:i w l ls 4 title '0.004',\
    # '../tr_char_p0002.csv' u 1:i w l ls 5 title '0.002'
# plot '../tr_char_100meV_HDS22_2000it_2e-5.csv' title columnhead(i) with lines lc rgb "white",\
#     '' u 1:i w l ls 1 title 'No dissipation',\
#     '../tr_char_100meV_HDS22_2000it_2e-5_tR1e-12.csv' u 1:i w l ls 2 title '{/Symbol t}_R = 1E-12 s',\
#     '../tr_char_100meV_HDS22_2000it_2e-5_tM1e-12.csv' u 1:i w l ls 3 title '{/Symbol t}_M = 1E-12 s',\
    # '../tr_char_100meV_HDS22_2000it_2e-5_tF1e-12.csv' u 1:i w l ls 4 title '{/Symbol t}_{/Symbol g} = 1E-12 s'
# plot '../tr_char_100meV_HDS22_2000it.csv' u 1:i w l ls 1
# plot '../tr_char_UDS3.csv' title columnhead(i) with lines lc rgb "white",\
#     '' u 1:i w l ls 1 title '150x150',\
#     '../tr_char_UDS3_200x200.csv' u 1:i w l ls 2 title '200x200'
# plot '../tr_char.csv' u 1:2 w l ls 1, '' u 1:3 w l ls 2,\
#     '' u 1:4 w l ls 3, '' u 1:5 w l ls 4, '' u 1:6 w l ls 5, '' u 1:7 w l ls 6
# plot '../tr_char_100meV_HDS22_2000it_2e-5.csv' title columnhead(i) with lines lc rgb "white",\
#     '' u 1:i w l ls 1 title '100meV',\
#     '../tr_char_10meV_HDS22_2000it_2e-5.csv' u 1:i w l ls 2 title '10meV',\
    # '../tr_char_0meV.csv' u 1:i w lp ls 3 title '0meV'
# plot '../tr_char_UDS2.csv' title columnhead(i) with lines lc rgb "white",\
#     '../tr_char_UDS2.csv' u 1:i w l ls 1 title 'UDS2',\
#     '../tr_char_UDS3.csv' u 1:i w l ls 2 title 'UDS3',\
#     '../tr_char_HDS22.csv' u 1:i w l ls 3 title 'HDS22'
