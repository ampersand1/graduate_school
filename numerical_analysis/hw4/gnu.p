!set terminal postscript enhanced
!set logscale x
!set xrange[0.001:0.1]

set title 'Comparison Between ||E|| and k for h=0.001'

set xlabel font ",12"
set ylabel font ",12"
set title font ",12"

set ylabel '||E||'
set xlabel 'k'

plot "1.0.001.dat" using 3:4 smooth unique title 'Jacobi',\
    "2.0.001.dat" using 3:4 smooth unique title 'G-S' , \
    "3.0.001.dat" using 3:4 smooth unique title 'SOR w-opt' , \
    "3.0.001-2.dat" using 3:4 smooth unique title 'SOR my w'

pause -1 "Hit any key to continue"
