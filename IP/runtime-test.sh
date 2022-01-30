#!/bin/bash

mkdir -p model

while getopts "frcdeanbhs" opt; do
	case $opt in
		f)	fixcol="_fixcol"
			f="-f" ;;
		r)	opts="_randobj" ;;
		c)	opts="_cuts" ;;
		h)	opts="_hints" ;;
		d)	fixcol="_fixdiag1" ;;
		e)	fixcol="_fixdiag2" ;;
		a)	fixcol="_fixantidiag" ;;
		n)	fixcol="_nofixcol"
			f="-n" ;;
		b)	fixcol="_backcirc"
			f="-n" ;;
		s)	fixcol="_nosym"
			f="-s" ;;
	esac
done
shift $(($OPTIND-1))

if [ -z $2 ]
then
	echo "Need order of the mutually orthogonal Latin squares, # of squares to find, and (optionally) # of tests to run"
	echo "Usage: $0 [-f] [-r] [-c] [-d] [-e] [-a] [-n] [-b] [-h] [-s] n k t"
	echo "-f fixes the first column of the second square"
	echo "-r enables a random objective function"
	echo "-c enables clique cuts via a custom callback function"
	echo "-d fixes the diagonal of the first square to contain zeros"
	echo "-e fixes the diagonal of the second square to contain zeros"
	echo "-a fixes the antidiagonal of the first square to contain (n-1)s"
	echo "-n disables fixing the first column in either square"
	echo "-b fixes the first square to be back-circulant"
	echo "-h uses variable hints"
	echo "-s uses no symmetry breaking"
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
make ${k}mols${fixcol}${opts}_$n

echo "time" > log/${k}mols${fixcol}${opts}_$n.$nt.times

for s in `shuf -i 0-2000000000 -n $nt`
do
	./${k}mols${fixcol}${opts}_$n $s | tee log/${k}mols${fixcol}${opts}_$n.$s.log
	t=`grep Explored log/${k}mols${fixcol}${opts}_$n.$s.log | cut -d' ' -f8`
	echo $t >> log/${k}mols${fixcol}${opts}_$n.$nt.times
	if grep "Optimal solution found" log/${k}mols${fixcol}${opts}_$n.$s.log -q
	then
		tail -n $((n*k+k-1)) log/${k}mols${fixcol}${opts}_$n.$s.log | python check.py $f $n $k
	fi
	numnodes=`grep Explored log/${k}mols${fixcol}${opts}_$n.$s.log | cut -d' ' -f2`
	echo "$numnodes/$t" | bc -l > log/${k}mols${fixcol}${opts}_$n.$nt.nodespeed
done

datamash max 1 min 1 mean 1 -H < log/${k}mols${fixcol}${opts}_$n.$nt.times
