
set title 'log h vs log error'
set logscale x
set logscale y
set ylabel 'log ||E||'
set xlabel 'log h'

set xlabel font ",12"
set ylabel font ",12"
set title font ",12"

f(x) = a*x+b
fit f(x) "hw2.dat" using (log($1)):(log($3)) via a, b


plot "hw2.dat" using 1:3 smooth unique title 'E_inf', "hw2.dat" using 1:4 smooth unique title 'E1', "hw2.dat" using 1:5 smooth unique title 'E2'

pause -1 "Hit any key to continue"




