#!/bin/bash

rm results.csv
touch results.csv

START=1
STOP=10000000

for (( k=10000000; k<=1000000000; k+=10000000 ));
do
	echo "acc,$START,$STOP,$k,$(./itg_acc.pgcc $START $STOP $k)">>results.csv;
	echo "ker,$START,$STOP,$k,$(./itg_ker.pgcc $START $STOP $k)">>results.csv;
	echo "omp,$START,$STOP,$k,$(./itg_omp.icc $START $STOP $k)">>results.csv;
done
