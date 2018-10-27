set terminal epslatex color size 14cm,14cm
set palette rgb 33,13,10
set output "Saurel1.tex"

set multiplot layout 2,2

set xrange [0:1]
set border lw 3
set size square
unset key
set xlabel '$x$'
set xtics 0.5
set tics nomirror
set ylabel rotate

set ylabel '$\rho$'
plot "./output/Saurel1.dat" u 1:2 w l lw 8


set ylabel '$u$'
plot "./output/Saurel1.dat" u 1:3 w l lw 8

set ylabel '$p$'
plot "./output/Saurel1.dat" u 1:4 w l lw 8

set ylabel '$e$'
plot "./output/Saurel1.dat" u 1:5 w l lw 8
