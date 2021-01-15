set term png
set output "ex_galerkin.png"
set size ratio -1
set xrange[0:1]
set yrange[0.2:-0.15]
set grid
set key top
set xlabel "x"
set ylabel "deflection"
p -(x**2-x)/2 lw 2 t "Exact sol.",\
"ex_galerkin.dat" lw 2 w lp t "Weak form + Galerkin."