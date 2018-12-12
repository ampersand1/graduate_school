!set terminal postscript enhanced
set logscale x
!set xrange[0.001:0.1]

set title 'Comparison Between ||E||/||E0|| and k for h=0.001'

set xlabel font ",12"
set ylabel font ",12"
set title font ",12"

set ylabel '||E||/||E0||'
set xlabel 'k'

plot "cg2.dat" using 1:2 smooth unique title 'cg',\
    "pcg2.dat" using 1:2 smooth unique title 'pcg'
pause -1 "Hit any key to continue"
