#!/bin/bash

while getopts "fbideanl" opt; do
	case $opt in
		f)	fixcol="_fixcol"
			f="-f" ;;
		d)	fixcol="_fixdiag"
			f="-f" ;;
		i)	model="_index" ;;
		b)	model="_bool" ;;
		e)	model="_extra" ;;
		a)	fixcol="_fixantidiag" ;;
		n)  fixcol="_nosym"
			f="-n" ;;
		l)  model="_linearelem" ;;
	esac
done
shift $(($OPTIND-1))

if [ -z $2 ]
then
	echo "Need order of the mutually orthogonal Latin squares and # of squares to find"
	echo "Usage: $0 [-f] [-d] [-a] [-i] [-b] [-e] [-n] [-l] n k"
	echo "-f fixes the first column of the second square"
	echo "-d additionally fixes both diagonals of the first square"
	echo "-a fixes the diagonal and anti-diagonal of the first square"
	echo "-i uses an encoding based on indexing constraints"
	echo "-b uses an encoding based on boolean constraints"
	echo "-e uses a boolean encoding with extra constraints"
	echo "-l uses a linear encoding based on element constraints"
	echo "-n uses no symmetry breaking"
	echo "n is the order of the Latin squares"
	echo "k is the number of MOLS to find (k=3 always uses indexing constraints)"
	exit
fi

n=$1
k=$2

mkdir -p log
make ${k}mols${fixcol}${model}_$n

echo "time" > log/${k}mols${fixcol}${model}_$n.time

./${k}mols${fixcol}${model}_$n $s | tee log/${k}mols${fixcol}${model}_$n.log
t=`grep elapsed log/${k}mols${fixcol}${model}_$n.log | cut -d' ' -f3`
echo $t >> log/${k}mols${fixcol}${model}_$n.time
if grep "Solution found" log/${k}mols${fixcol}${model}_$n.log -q
then
	python check.py $n $k $f < log/${k}mols${fixcol}${model}_$n.log
fi
