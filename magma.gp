set terminal post eps color enhanced size 11cm,6cm
set out "magma.eps"

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

plot [0:190][0:1000] 'test_irred.dat' using 2:($6+$7) with line lt 1 lw 3 title "irred",\
'test_ext.dat' using 2:($5+$6) with line lt 2 lw 3 title "P R",\
'test_compositum.dat' using 2:($5+$6+$7) with line lt 3 lw 3 title "P Q",\
'test_relative.dat' using 2:($5+$6) with line lt 4 lw 3 title "ext"
#'bench.dat' index 1 using 1:8 with line lt 2 lw 3 title "mulmod",\
#'bench.dat' index 1 using 1:9 with line lt 3 lw 3 title "tmulmod",\
#'bench.dat' index 1 using 1:4 with line lt 4 lw 3 title "Embed",\
#'bench.dat' index 1 using 1:($5-$6) with line lt 5 lw 3 title "D2M"

unset ylabel
set size 0.4,1
set origin 0.55,0
plot [0:190] 'test_relative.dat' using 2:($7/1000) with line lt 1 lw 3 title "embedding",\
'test_relative.dat' using 2:($8/1000) with line lt 2 lw 3 title "inverse isomorphism"


unset multiplot