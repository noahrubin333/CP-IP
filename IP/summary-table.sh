#!/bin/bash

while getopts "frcdeanbhs" opt; do
	case $opt in
		f) fixcol="_fixcol" ;;
		r) opts="_randobj" ;;
		c) opts="_cuts" ;;
		d) fixcol="_fixdiag1" ;;
		e) fixcol="_fixdiag2" ;;
		a) fixcol="_fixantidiag" ;;
		n) fixcol="_nofixcol" ;;
		b) fixcol="_backcirc" ;;
		h) opts="_hints" ;;
		s) fixcol="_nosym" ;;
	esac
done
shift $(($OPTIND-1))

if [ -z $1 ]
then
	echo "Need # of squares searched for and (optionally) the # of tests that were run"
	exit
fi

k=$1

if [ -z $2 ]
then
	c=1
else
	c=$2
fi

printf "order\tavg\tmin\tmax\tnodes/sec\n"

for n in `seq 3 12`
do
	if [ ! -f log/${k}mols${fixcol}${opts}_$n.$c.times ]
	then
		continue
	fi
	printf "$n\t"
	t=`datamash mean 1 --header-in < log/${k}mols${fixcol}${opts}_$n.$c.times`
	printf "%.1f\t" $t
	t=`datamash min 1 --header-in < log/${k}mols${fixcol}${opts}_$n.$c.times`
	printf "%.1f\t" $t
	t=`datamash max 1 --header-in < log/${k}mols${fixcol}${opts}_$n.$c.times`
	printf "%.1f\t" $t
	if [ -f log/${k}mols${fixcol}${opts}_$n.$c.nodespeed ]
	then
		nodespeed=`cat log/${k}mols${fixcol}${opts}_$n.$c.nodespeed`
		printf "%.2f" $nodespeed
	fi
	printf "\n"
done

