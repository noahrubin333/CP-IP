#!/bin/bash

while getopts "fntcd" opt; do
	case $opt in
		f)	cuts="_fixcol"
			f="-f" ;;
		n)	cuts="_nofixcol"
			f="-n" ;;
		t)	cuts="_fixcolthird"
			f="-n" ;;
		c)	cuts="_cycletype"
			f="-n" ;;
		d)	cuts="_domainred"
			f="-n" ;;
	esac
done
shift $(($OPTIND-1))

if [ -z $2 ]
then
	echo "Need order of the mutually orthogonal Latin squares and # of squares to find"
	echo "Usage: $0 [-f] [-n] [-t] [-d] [-c] n k"
	echo "-n performs no fixing of the first column"
	echo "-f performs fixing in the first column"
	echo "-t fixes the first column of the third square to appear in sorted order"
	echo "-d uses domain reduction on the first column of the second square (and implies -t)"
	echo "-c ensures the first column of the second square is one of a fixed set of cycle type representatives (and implies -t)"
	echo "n is the order of the Latin squares"
	echo "k is the number of Latin squares to find"
	exit
fi

n=$1
k=$2
#t=$3

mkdir -p log
make ${k}mols${cuts}_$n

./${k}mols${cuts}_$n | tee log/${k}mols${cuts}_$n.log
if grep "Solution found" log/${k}mols${cuts}_$n.log -q
then
	command="tail -n $((k*(n+2)+1)) log/${k}mols${cuts}_$n.log | python check.py $n $k $f"
	#echo $command
	eval $command
fi
