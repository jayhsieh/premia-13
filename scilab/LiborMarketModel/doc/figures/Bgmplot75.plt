set terminal postscript enhanced color
set output "SwaptionPayeur75.eps";
set title "Payer Swaption, K=0.075,TS=4,r=0.1,{/Symbol t=1}"
set key right bottom Left width 1 box
set xr[4:34]
set xlabel "Nb Steps"
plot "resultat1_0_0_4_3_2_075_4_1.data" u 1:3 title "-Bushy-Rotate" w l ,"resultat1_1_0_4_3_2_075_4_1.data" u 1:3 title "+Bushy-Rotate" w l ,"resultat1_0_1_4_3_2_075_4_1.data" u 1:3 title "-Bushy+Rotate" w l,"resultat1_1_1_4_3_2_075_4_1.data" u 1:3 title "+Bushy+Rotate" w l
set terminal windows

set terminal postscript enhanced color
set output "SwaptionReceveur75.eps";
set title "Payer Swaption, K=0.075,TS=4,r=0.1,{/Symbol t=1}"
set key right bottom Left width 1 box
set xr[4:34]
set xlabel "Nb Steps"
plot "resultat_moin1_0_0_4_3_2_075_4_1.data" u 1:3 title "-Bushy-Rotate" w l ,"resultat_moin1_1_0_4_3_2_075_4_1.data" u 1:3 title "+Bushy-Rotate" w l ,"resultat_moin1_0_1_4_3_2_075_4_1.data" u 1:3 title "-Bushy+Rotate" w l,"resultat_moin1_1_1_4_3_2_075_4_1.data" u 1:3 title "+Bushy+Rotate" w l
set terminal windows
