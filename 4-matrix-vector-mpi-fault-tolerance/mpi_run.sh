#!/bin/bash -eu
mpicc mpi.c -o mpi
mpiexec -np $NP --with-ft ulfm mpi $SIZE
