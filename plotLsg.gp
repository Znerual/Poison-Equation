set term pdfcairo
set output "lsg.pdf"

set key font ",6"
set xtics font ",6"
set ytics font ",6"

set xlabel "x_i"
set ylabel "u_i"
set title "Potential n = 10"
plot "lsgn10.dat" using 1:2 title "Exakt" w lp, "lsgn10.dat" using 1 : 3 title "Gauﬂ" w lp, "lsgn10.dat" using 1 : 4 title "Thomas" w lp, "lsgn10.dat" using 1 : 5 title "DGESV" w lp
#plot "lsgn10.dat" using 1:2 title "Exakt" w lp, "lsgn10.dat" using 1 : 3 title "Gauﬂ" w lp, "lsgn10.dat" using 1 : 4 title "Thomas" w lp, "lsgn10.dat" using 1 : 5 title "DGESV" w lp, "lsgn10.dat" using 1 : 6 title "DPTSV" w lp
set title "Potential n = 100"
plot "lsgn100.dat" using 1:2 title "Exakt" w lp, "lsgn100.dat" using 1 : 3 title "Gauﬂ" w lp, "lsgn100.dat" using 1 : 4 title "Thomas" w lp, "lsgn100.dat" using 1 : 5 title "DGESV" w lp
#plot "lsgn100.dat" using 1:2 title "Exakt" w lp, "lsgn100.dat" using 1 : 3 title "Gauﬂ" w lp, "lsgn100.dat" using 1 : 4 title "Thomas" w lp, "lsgn100.dat" using 1 : 5 title "DGESV" w lp, "lsgn100.dat" using 1 : 6 title "DPTSV" w lp

unset out
unset term