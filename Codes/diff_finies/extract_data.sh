#!/bin/bash

rm results.csv
touch results.csv

ITER=1000
for (( k=1000; k<=6000; k+=200 ));
do
	echo "acc,$k,$ITER,$(./diff_finies_acc.pgcc $k $ITER)">>results.csv;
	echo "ker,$k,$ITER,$(./diff_finies_ker.pgcc $k $ITER)">>results.csv;
	echo "omp,$k,$ITER,$(./diff_finies_omp.icc $k $ITER)">>results.csv;
done
