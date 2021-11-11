#!/bin/bash

alias progess=echo "$(cat results.csv | wc -l)"

rm results.csv
touch results.csv

for (( i=400; i<=7000; i+=200 ));
do
	echo "omp,$i,$(./cblas_omp.icc $i)">>results.csv;
	echo "ker,$i,$(./cblas_ker.pgcc $i)">>results.csv;
	echo "acc,$i,$(./cblas_acc.pgcc $i)">>results.csv;
done;	
