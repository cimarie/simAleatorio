set terminal png medium font Arial
set output "Nequilibrio_beta=0.5.png"
set title "Equilibria stability"
set xlabel "migration rate"
set ylabel "initial fraction of altruists"
set autoscale xy
set linestyle 1 lt 9 lw 1
set key box linestyle 1
set key center
plot "1eq_beta=0.5" using 1:2 with points pt 6 ps 2 lt 7 title "stable", "2eq_beta=0.5" using 1:2 with points pt 7 ps 2 lt 7 title "unstable" 
