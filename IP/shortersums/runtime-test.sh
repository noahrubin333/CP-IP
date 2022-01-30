#!/bin/bash

while getopts "f" opt; do
	case $opt in
		f)	fixcol="_fixcol"
			f="-f" ;;
	esac
done
shift $(($OPTIND-1))

if [ -z $2 ]
then
	echo "Need order of the mutually orthogonal Latin squares, # of squares to find, and (optionally) # of tests to run"
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
make ${k}mols${fixcol}_$n

echo "time" > log/${k}mols${fixcol}_$n.$nt.times

for s in `shuf -i 0-2000000000 -n $nt`
do
	./${k}mols${fixcol}_$n $s | tee log/${k}mols${fixcol}_$n.$s.log
	t=`grep Explored log/${k}mols${fixcol}_$n.$s.log | cut -d' ' -f8`
	echo $t >> log/${k}mols${fixcol}_$n.$nt.times
	if grep "Optimal solution found" log/${k}mols${fixcol}_$n.$s.log -q
	then
		tail -n $((n*k+k-1)) log/${k}mols${fixcol}_$n.$s.log | python check.py $f $n $k
	fi
done

datamash max 1 min 1 mean 1 -H < log/${k}mols${fixcol}_$n.$nt.times
