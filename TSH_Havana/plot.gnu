p "salida" u 1:2 w l t "poblacion-27","salida" u 1:3 w l  t "poblacion-26","salida" u 1:4 w lp t "acoplamiento real","salida" u 1:5 w lp t "acoplamiento imag","salida" u 1:6 w lp t "acoplamiento total"


set grid xtics nomxtics ytics nomytics noztics nomztics \
 nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   lt 0 linewidth 0.500,  lt 0 linewidth 0.500

set xlabel "t (ps) "
set ylabel "Porciento %"
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title "" center
set key out right top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 4 spacing 1.2 width 0 height 0
set key maxcolumns 0 maxrows 0
set terminal postscript eps enhanced color colortext
set output "Coupling.eps"
rep

