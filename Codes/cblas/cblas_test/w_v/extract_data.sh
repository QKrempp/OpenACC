#!/bin/bash

rm results.csv
touch results.csv
for i in *.pgcc;
do
	echo $i;
	GANG=$(echo $i | grep -oP "[0-9]{1,3}(?=g)");
	WORKER=$(echo $i | grep -oP "[0-9]{1,3}(?=w)");
	VECTOR=$(echo $i | grep -oP "[0-9]{1,3}(?=v)");
	for j in {1000..10000..1000};
	do
		echo "$GANG,$WORKER,$VECTOR,$j,$(./$i $j)">>results.csv;
	done
done
