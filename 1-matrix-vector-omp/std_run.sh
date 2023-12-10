#!/bin/bash -eu
xlc -qsmp=omp -o omp omp.c
for i in {0..6}
do
	echo ===$((125*2**$i))===
	for j in {0..6}
	do
		OMP_NUM_THREADS=$((2**$j)) ./omp $((125*2**$i))
	done
done
