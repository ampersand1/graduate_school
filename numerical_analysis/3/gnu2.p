! Solution graphing

set title 'Comparison Between Exact Solution and Approximation to e^{x}-x+(1-e) using h= 0.001'

set ylabel 'y'
set xlabel 'x'

set xlabel font ",12"
set ylabel font ",12"
set title font ",12"


plot exp(x)-x+(1-2.718281828459) title 'real', "method1-2.dat" using 1:2 smooth unique title 'method 1', "method2-2.dat" using 1:2 smooth unique title 'method 2',"method3-2.dat" using 1:2 smooth unique title 'method 3'


pause -1 "Hit any key to continue"
