#!/bin/bash

f="-p"

while getopts "unftcd" opt; do
	case $opt in
		u)	cuts="_cuts" ;;
		n)	cuts="_nofixcol"
			f="-n" ;;
		f)	cuts="_fixcol"
			f="-f" ;;
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
	echo "Need order of the mutually orthogonal Latin squares, the number of Latin squares to find, and optionally the # of transversals to fix using CP solver (in the 2MOLS model)"
	echo "Usage: $0 [-u] [-n] [-f] [-t] [-d] [-c] n k t"
	echo "-u enables clique cuts via a custom callback function (ignored in the 3MOLS model)"
	echo "-n performs no fixing of the first column"
	echo "-f fixes the first column of the second square"
	echo "-t fixes the first column of the third square to appear in sorted order"
	echo "-d uses domain reduction on the first column of the second square (and implies -t)"
	echo "-c ensures the first column of the second square is one of a fixed set of cycle type representatives (and implies -t)"
	echo "n is the order of the Latin squares"
	echo "k is the number of Latin squares to find"
	echo "t is the number of transversals to fix (only needed in the 2MOLS model)"
	exit
fi

n=$1
k=$2
t=$3

if [ "$k" == "3" ] || [ "$k" == "2.5" ]
then
	mkdir -p log
	make ${k}mols${cuts}_$n

	./${k}mols${cuts}_$n | tee log/${k}mols${cuts}_$n.log
	if [ "$k" == "3" ] && grep "Found solution" log/${k}mols${cuts}_$n.log -q
	then
		tail -n $((3*(n+1)+1)) log/${k}mols${cuts}_$n.log | head -n $((3*n+2)) | python check.py $f $n $k
	fi
else
	mkdir -p log
	make ${k}mols$cuts

	command="./${k}mols$cuts $n $t | tee log/${k}mols${cuts}_${n}_$t.log"
	echo $command
	eval $command
	if grep "Optimal solution found" log/${k}mols${cuts}_${n}_$t.log -q
	then
		tail -n $((2*n+1)) log/${k}mols${cuts}_${n}_$t.log | python check.py $f $n $k
	fi
fi
