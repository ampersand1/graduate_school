
!set terminal postscript enhanced
set title 'Comparison Between log ||E|| and log(h) for Method 3'



f(x) = a*x+b

! First is for infinity, second is for l2.
fit [log(0.005):log(0.1)] f(x) "method3.dat" using (log($1)):(log($3)) via a, b
fit [log(0.005):log(0.1)] f(x) "method3.dat" using (log($1)):(log($4)) via a, b

set xrange[0.001:0.1]
set logscale x
set logscale y

set ylabel 'log ||E||'
set xlabel 'log(h)'

set xlabel font ",12"
set ylabel font ",12"
set title font ",12"


plot "method3.dat" using 1:3 smooth unique title '||E_{infty}||', "method3.dat" using 1:4 smooth unique title '||E_{L^2}||'


pause -1 "Hit any key to continue"
