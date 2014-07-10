#!/bin/bash
#Usage : ./test.sh test_number fixed_parameters number_of_tests

N=1
SAGE=${SAGE:-sage}
GNUPLOT=${GNUPLOT:-gnuplot}
CMD=batch.py


if [[ "$1" == 1 ]];then
    while [ -f "data_test1_p_fixed${N}.txt"  ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test1_p_fixed${N}"

    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 > ${NAMEFILE}.txt

    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"extension degree\"
    set ylabel \"time(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'time_m' with lines,\
     \"${NAMEFILE}.txt\" using 1:3 title 'time_E' with lines,\
     \"${NAMEFILE}.txt\" using 1:4 title 'time_ordm' with lines,\
     \"${NAMEFILE}.txt\" using 1:5 title 'time_period' with lines"| $GNUPLOT

elif [[ "$1" == 2 ]];then
    while [ -f "data_test2_p_fixed${N}.txt"  ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test2_p_fixed${N}.txt"

    echo 'TODO'


elif [[ "$1" == 3 ]];then
    while [ -f "data_test3_n_fixed${N}.txt"  ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test3_n_fixed${N}"

    echo 'Generating the data file...'
    SAGE $CMD $1 $2 $3 > ${NAMEFILE}.txt

    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"characteristic\"
    set ylabel \"time(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'time_m' with lines,\
     \"${NAMEFILE}.txt\" using 1:3 title 'time_E' with lines,\
     \"${NAMEFILE}.txt\" using 1:4 title 'time_ordm' with lines,\
     \"${NAMEFILE}.txt\" using 1:5 title 'time_period' with lines"| $GNUPLOT


elif [[ "$1" == 4 ]];then
    while [ -f "data_test4_n_fixed${N}.txt"  ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test4_n_fixed${N}"

    echo 'TODO'

fi
