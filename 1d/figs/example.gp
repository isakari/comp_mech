set term png
set output "example.png"
set size ratio -1
set xrange[0:1]
set yrange[0.2:-0.15]
set grid
set key top
set xlabel "x"
set ylabel "deflection"
p -(x**2-x)/2 lw 2 t "Exact sol.",\
  1/(2*pi)*sin(pi*x) lw 2 t "WRM with v=1",\
  4/(pi**3)*sin(pi*x) lw 2 t "WRM with the Galerkin method."