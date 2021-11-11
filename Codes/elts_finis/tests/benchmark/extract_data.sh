#!/bin/bash
module load gmsh
rm results.csv
rm square.geo

for i in *.geo
do
	echo $i
	MAILLE=$(echo $i|grep -oP "(?<=square)[0-9]{3}");
	gmsh $i -2;
	mv $(echo "square$MAILLE.msh") square.msh;
	echo "icc,$MAILLE,$(../myfem.icc)">>results.csv;
	echo "pgcc,$MAILLE,$(../myfem.pgcc)">>results.csv;
done;
