#!/bin/bash

while getopts "fbideanl" opt; do
	case $opt in
		f) fixcol="_fixcol" ;;
		d) fixcol="_fixdiag" ;;
		i) model="_index" ;;
		b) model="_bool" ;;
		e) model="_extra" ;;
		a) fixcol="_fixantidiag" ;;
		n) fixcol="_nosym" ;;
		l) model="_linearelem" ;;
	esac
done
shift $(($OPTIND-1))

if [ -z $1 ]
then
	echo "Need # of squares searched for"
	echo "Usage: $0 [-f] [-d] [-a] [-i] [-b] [-e] [-n] k"
	echo "-f fixes the first column of the second square"
	echo "-d fixes the first column of both squares and a diagonal of each"
	echo "-a fixes the diagonal and anti-diagonal of the first square"
	echo "-i encodes orthogonality as an indexing constraint"
	echo "-b encodes orthogonality as boolean constraints"
	echo "-e uses boolean encoding with some extra constraints"
	echo "-l uses a linear encoding based on element constraints"
	echo "-n uses no symmetry breaking"
	echo "k is the number of MOLS searched for"
	exit
fi

k=$1

printf "order\tavg time\n"

for n in `seq 3 10`
do
	if [ ! -f log/${k}mols${fixcol}${model}_$n.time ] || [ $(wc -l log/${k}mols${fixcol}${model}_$n.time | cut -d' ' -f1) -eq 1 ]
	then
		continue
	fi
	printf "$n\t"
	t=`datamash mean 1 --header-in < log/${k}mols${fixcol}${model}_$n.time`
	printf "%.1f\n" `echo "$t/1000" | bc -l`
done

