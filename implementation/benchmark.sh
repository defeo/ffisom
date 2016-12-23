#! /bin/bash
SAGE=${SAGE:-sage}
fname=$1.dat
echo -ne "\033]0;$0 $*\007"
echo "# $0 $*" >> $fname
echo "# `cat /proc/cpuinfo | grep "model" | tail -1`" >> $fname
echo "# `$SAGE --version`" >> $fname
echo "# ffisom revision `git rev-parse HEAD`" >> $fname
echo "# Magma `$SAGE -c 'print(magma.version()[1])' 2> /dev/null`" >> $fname
command="import os; sys.path.append(os.getcwd()); from benchmark import *; benchmark(pbound=[${2:-3}, ${3:-Infinity}], nbound=[${4:-3}, ${5:-Infinity}], prime=${6:-2}, loops=${7:-5}, tloop=${8:-Infinity}, verbose=${9:-True}, write=True, overwrite=False, fname=os.getcwd()+'/$fname', ${10})"
echo "# command=\"$command\"" >> $fname
LD_LIBRARY_PATH=. $SAGE -c $command
