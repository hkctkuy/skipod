#!/bin/bash -eu
mpicc mpi.c -o mpi
mpiexec -np $NP mpi $SIZE
