#!/bin/bash

while getopts "f" opt; do
	case $opt in
		f) fixcol="_fixcol" ;;
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

printf "order\tavg\tmin\tmax\n"

for n in `seq 3 10`
do
	if [ ! -f log/${k}mols${fixcol}_$n.$c.times ]
	then
		continue
	fi
	printf "$n\t"
	t=`datamash mean 1 --header-in < log/${k}mols${fixcol}_$n.$c.times`
	printf "%.1f\t" $t
	t=`datamash min 1 --header-in < log/${k}mols${fixcol}_$n.$c.times`
	printf "%.1f\t" $t
	t=`datamash max 1 --header-in < log/${k}mols${fixcol}_$n.$c.times`
	printf "%.1f\n" $t
done

