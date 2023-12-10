#!/bin/bash -eu
mpicc scan.c -o scan
mpiexec -np 16 --oversubscribe ./scan
