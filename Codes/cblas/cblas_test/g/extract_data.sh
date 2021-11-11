#!/bin/bash

rm results.csv
touch results.csv
for i in *.pgcc;
do
	echo $i;
	GANG=$(echo $i | grep -oP "[0-9]{1,5}(?=g)");
	WORKER=$(echo $i | grep -oP "[0-9]{1,4}(?=w)");
	VECTOR=$(echo $i | grep -oP "[0-9]{1,4}(?=v)");
	for (( k=0; k<=3; k++ ));
	do
		if [ $GANG -lt 128 -a $GANG -ne 0 ];
		then
			for (( j=512; j<=8192; j=$j*2 ));
			do
				echo "$GANG,$WORKER,$VECTOR,$j,$(./$i $j)">>results.csv;
			done
		else
			for (( j=512; j<=16384; j=$j*2 ));
			do
				echo "$GANG,$WORKER,$VECTOR,$j,$(./$i $j)">>results.csv;
			done
		fi
	done
done
