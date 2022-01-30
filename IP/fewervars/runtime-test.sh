#!/bin/bash

while getopts "fne" opt; do
	case $opt in
		f)	fixcol="_fixcol"
			f="-f" ;;
		n)	vars="_nonewvars" ;;
		e)	vars="_extravars" ;;
	esac
done
shift $(($OPTIND-1))

if [ -z $2 ]
then
	echo "Need order of the mutually orthogonal Latin squares, # of squares to find, and (optionally) # of tests to run"
	echo "Usage: $0 [-f] [-n] [-e] n k t"
	echo "-f fixes the first column of the second square"
	echo "-n uses a 2MOLS(n) model with only 2n^3 variables"
	echo "-e uses a 2MOLS(n) model with an extra n^3 variables"
	echo "n is the order of the Latin squares"
	echo "k is the number of MOLS to find"
	echo "t is the number of tests to run"
	exit
fi

n=$1
k=$2

if [ -z $3 ]
then
	nt=1
else
	nt=$3
fi

mkdir -p log
make ${k}mols${fixcol}${vars}_$n

echo "time" > log/${k}mols${fixcol}${vars}_$n.$nt.times

for s in `shuf -i 0-2000000000 -n $nt`
do
	./${k}mols${fixcol}${vars}_$n $s | tee log/${k}mols${fixcol}${vars}_$n.$s.log
	t=`grep Explored log/${k}mols${fixcol}${vars}_$n.$s.log | cut -d' ' -f8`
	echo $t >> log/${k}mols${fixcol}${vars}_$n.$nt.times
	if grep "Optimal solution found" log/${k}mols${fixcol}${vars}_$n.$s.log -q
	then
		tail -n $((n*k+k-1)) log/${k}mols${fixcol}${vars}_$n.$s.log | python check.py $f $n $k
	fi
done

datamash max 1 min 1 mean 1 -H < log/${k}mols${fixcol}${vars}_$n.$nt.times
