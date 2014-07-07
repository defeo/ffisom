N=1

while [ -f "data${N}.txt"  ]
do
    N=$((N+1))
done

NAMEFILE="data${N}.txt"

echo 'Generating the data file...'
sage print.sage $1 $2 $3 > $NAMEFILE

echo 'Generating the graph...'

echo "plot \"$NAMEFILE\" using 1:2 title 'time_m' with lines, \"$NAMEFILE\"\
using 1:3 title 'time_E' with lines, \"$NAMEFILE\" using 1:4 title\
'time_ordm' with lines, \"$NAMEFILE\" using 1:5 title 'time_period' with\
 lines" | gnuplot
