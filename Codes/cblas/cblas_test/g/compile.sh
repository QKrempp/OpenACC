#!/bin/bash

for i in *.c;
do
	RAD=$(echo $i| grep -oP ".*(?=.c)");
	pgcc -fast -acc "$RAD".c -o "$RAD".pgcc;
done

