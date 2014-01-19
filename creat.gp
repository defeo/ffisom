set terminal post eps color enhanced size 11cm,6cm
set out "creat.eps"

set key left top
set rmargin 0
set lmargin 0
#set size square

set multiplot 

set size 0.425,1
set origin 0.07,0
#set ylabel "seconds"
set xlabel "m"
set xtics 0, 30

plot [0:190][0:0.08] 'bench.dat' index 1 using 1:2 with line lt 1 lw 3 title "R",\
'bench.dat' index 1 using 1:8 with line lt 2 lw 3 title "mulmod",\
'bench.dat' index 1 using 1:9 with line lt 3 lw 3 title "tmulmod",\
'bench.dat' index 1 using 1:4 with line lt 4 lw 3 title "Embed",\
'bench.dat' index 1 using 1:($5-$6) with line lt 5 lw 3 title "D2M",


unset ylabel
set size 0.425,1
set origin 0.55,0
plot [0:190][0:3] 'bench.dat' index 1 using 1:2 with line lt 1 lw 3 title "R",\
'bench.dat' index 1 using 1:10 with line lt 8 lw 3 title "Phi1",\
'bench.dat' index 1 using 1:11 with line lt 9 lw 3 title "Phi2"


unset multiplot